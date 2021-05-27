#!/usr/bin/env python3
# Roger Volden and Chris Vollmers

import sys
import os
import argparse
import mappy as mm
from tqdm import tqdm
import multiprocessing as mp
import editdistance as ld
from glob import glob
import gzip
import shutil

VERSION = 'v2.2.3'

def parse_args():
    '''Parses arguments.'''
    parser = argparse.ArgumentParser(description='Reorients/demuxes/trims consensus reads.',
                                     add_help=True,
                                     prefix_chars='-')
    parser.add_argument('--input_fasta_file', '-i', type=str, action='store',
                        help='Fasta file with consensus called R2C2 reads')
    parser.add_argument('--output_path', '-o', type=str, action='store', default=os.getcwd(),
                        help='''Directory where all the files will end up.
                                Defaults to your current directory.''')
    parser.add_argument('--adapter_file', '-a', type=str, action='store',
                        help='Fasta file with adapter (3 and 5 prime) sequences')
    parser.add_argument('--index_file', '-x', type=str, action='store',
                        help='Fasta file with oligo dT indexes')
    parser.add_argument('--config', '-c', type=str, action='store', default='',
                        help='If you want to use a config file to specify paths to\
                              programs, specify them here. Use for poa, racon, water,\
                              blat, and minimap2 if they are not in your path.')
    parser.add_argument('--undirectional', '-u', action='store_true',
                        help='''By default, your cDNA molecules are assumed to be
                                directional with two sequences named "3Prime_adapter"
                                and "5Prime_adapter" expected in your adapter_file in
                                fasta format. If you add this flag your cDNA molecules
                                are expected to be undirectional and only one sequence
                                named "Adapter" should be in your adapter_file in fasta
                                format''')
    parser.add_argument('--trim', '-t', action='store_true',
                        help='Use this flag to trim the adapters off the ends of \
                              your sequences.')
    parser.add_argument('--barcoded', '-b', action='store_true', default=False,
                        help='Use if postprocessing 10x reads. Produces a separate \
                              file with 10x barcode sequences')
    parser.add_argument('--threads', '-n', type=int, default=1,
                        help='Number of threads to use during multiprocessing. Defaults to 1.')
    parser.add_argument('--groupSize', '-g', type=int, default=1000,
                        help='Number of reads processed by each thread in each iteration. Defaults to 1000.')
    parser.add_argument('--blatThreads', '-bt', action='store_true', default=False,
                        help='''Use to chunk blat across the number of threads instead of by groupSize (faster).''')
    parser.add_argument('--compress_output', '-co', action='store_true', default=False,
                        help='Use to compress (gzip) both the consensus fasta and subread fastq output files.')
    parser.add_argument('--version', '-v', action='version', version=VERSION, help='Prints the C3POa version.')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()

def configReader(path, configIn):
    progs = {}
    with open(configIn) as f:
        for line in f:
            if line.startswith('#') or not line.rstrip().split():
                continue
            line = line.rstrip().split('\t')
            progs[line[0]] = line[1]
    possible = set(['racon', 'blat'])
    inConfig = set()
    for key in progs.keys():
        inConfig.add(key)
    # check for missing programs
    # if missing, default to path
    for missing in possible-inConfig:
        path = missing
        progs[missing] = path
        sys.stderr.write('Using ' + str(missing)
                         + ' from your path, not the config file.\n')
    return progs

def get_file_len(inFile):
    '''Figure out how many reads for best chunk size for parallelization'''
    count = 0
    for _ in mm.fastx_read(inFile, read_comment=False):
        count += 1
    return count

def cat_files(path, pattern, output, pos, compress):
    '''Use glob to get around bash argument list limitations'''
    if compress:
        output += '.gz'
        final_fh = gzip.open(output, 'wb+')
    else:
        final_fh = open(output, 'w+')
    for f in tqdm(glob(path + pattern), position=pos):
        with open(f) as fh:
            for line in fh:
                if compress:
                    line = line.encode()
                final_fh.write(line)
    final_fh.close()

def remove_files(path, pattern):
    '''Use glob to get around bash argument list limitations'''
    for d in tqdm(glob(path + pattern), desc='Removing files'):
        shutil.rmtree(d)

def process(args, reads, blat, iteration, idx_to_seq, seq_to_idx):
    tmp_dir = args.output_path + 'post_tmp_' + str(iteration) + '/'
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)
    tmp_fa = tmp_dir + 'tmp_for_blat.fasta'
    tmp_fa_fh = open(tmp_fa, 'w+')
    for header, seq in reads.items():
        print('>' + header, file=tmp_fa_fh)
        print(seq, file=tmp_fa_fh)
    tmp_fa_fh.close()

    run_blat(tmp_dir, tmp_fa, args.adapter_file, blat)
    os.remove(tmp_fa)
    adapter_dict = parse_blat(tmp_dir, reads)
    write_fasta_file(args, tmp_dir, adapter_dict, reads, seq_to_idx, idx_to_seq)

def chunk_process(num_reads, args, blat):
    '''Split the input fasta into chunks and process'''
    if args.blatThreads:
        chunk_size = (num_reads // args.threads) + 1
    else:
        chunk_size = args.groupSize
    if chunk_size > num_reads:
        chunk_size = num_reads

    if args.index_file:
        idx_to_seq, seq_to_idx = read_fasta(args.index_file, True)
    else:
        idx_to_seq, seq_to_idx = {}, {}

    pool = mp.Pool(args.threads)
    pbar = tqdm(total=num_reads // chunk_size + 1, desc='Aligning with BLAT and processing')
    iteration, current_num, tmp_reads, target = 1, 0, {}, chunk_size
    for read in mm.fastx_read(args.input_fasta_file, read_comment=False):
        tmp_reads[read[0]] = read[1]
        current_num += 1
        if current_num == target:
            pool.apply_async(
                process,
                args=(args, tmp_reads, blat, iteration, idx_to_seq, seq_to_idx),
                callback=lambda _: pbar.update(1)
            )
            iteration += 1
            target = chunk_size * iteration
            if target >= num_reads:
                target = num_reads
            tmp_reads = {}
    pool.close()
    pool.join()
    pbar.close()

    flc = 'R2C2_full_length_consensus_reads.fasta'
    flc_left = 'R2C2_full_length_consensus_reads_left_splint.fasta'
    flc_right = 'R2C2_full_length_consensus_reads_right_splint.fasta'
    pool = mp.Pool(args.threads)
    print('Catting files', file=sys.stderr)
    if idx_to_seq:
        idx_to_seq['no_index_found'] = ''
        for idx in idx_to_seq.keys():
            idx += '/'
            if not os.path.isdir(args.output_path + idx):
                os.mkdir(args.output_path + idx)
            pattern = 'post_tmp*/' + idx
            pool.apply_async(cat_files, args=(
                    args.output_path,
                    pattern + flc,
                    args.output_path + idx + flc,
                    0, args.compress_output))
            pool.apply_async(cat_files, args=(
                    args.output_path,
                    pattern + flc_left,
                    args.output_path + idx + flc_left,
                    1, args.compress_output))
            pool.apply_async(cat_files, args=(
                    args.output_path,
                    pattern + flc_right,
                    args.output_path + idx + flc_right,
                    2, args.compress_output))
        mux_tsvs = 'post_tmp*/R2C2_oligodT_multiplexing.tsv'
        mux_tsv_final = args.output_path + 'R2C2_oligodT_multiplexing.tsv'
        pool.apply_async(cat_files, args=(args.output_path, mux_tsvs, mux_tsv_final, 3, False))
    else:
        pattern = 'post_tmp*/'
        pool.apply_async(cat_files, args=(
                args.output_path,
                pattern + flc,
                args.output_path + flc,
                0, args.compress_output))
        pool.apply_async(cat_files, args=(
                args.output_path,
                pattern + flc_left,
                args.output_path + flc_left,
                1, args.compress_output))
        pool.apply_async(cat_files, args=(
                args.output_path,
                pattern + flc_right,
                args.output_path + flc_right,
                2, args.compress_output))
        if args.barcoded:
            flc_bc = pattern + 'R2C2_full_length_consensus_reads_10X_sequences.fasta'
            flc_bc_final = args.output_path + 'R2C2_full_length_consensus_reads_10X_sequences.fasta'
            pool.apply_async(cat_files, args=(args.output_path, flc_bc, flc_bc_final, 3, args.compress_output))
    pool.close()
    pool.join()
    remove_files(args.output_path, 'post_tmp*')

def read_fasta(inFile, indexes):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readDict, index_dict = {}, {}
    for read in mm.fastx_read(inFile, read_comment=False):
        readDict[read[0]] = read[1]
        if indexes:
            index_dict[read[1]] = read[0]
    if indexes:
        return readDict, index_dict
    return readDict

def run_blat(path, infile, adapter_fasta, blat):
    align_psl = path + 'adapter_to_consensus_alignment.psl'
    if not os.path.exists(align_psl) or os.stat(align_psl).st_size == 0:
        os.system('{blat} -noHead -stepSize=1 -tileSize=6 -t=DNA -q=DNA -minScore=10 \
                  -minIdentity=10 -minMatch=1 -oneOff=1 {adapters} {reads} {psl} >{blat_msgs}'
                  .format(blat=blat, adapters=adapter_fasta, reads=infile, psl=align_psl, blat_msgs=path + 'blat_msgs.log'))
    else:
        print('Reading existing psl file', file=sys.stderr)

def parse_blat(path, reads):
    adapter_dict, iterator = {}, 0

    for name, sequence in reads.items():
        adapter_dict[name] = {}
        adapter_dict[name]['+'] = []
        adapter_dict[name]['-'] = []
        adapter_dict[name]['+'].append(('-', 1, 0))
        adapter_dict[name]['-'].append(('-', 1, len(sequence)))

    with open(path + 'adapter_to_consensus_alignment.psl') as f:
        for line in f:
            a = line.strip().split('\t')
            read_name, adapter, strand = a[9], a[13], a[8]
            if int(a[5]) < 50 and float(a[0]) > 10:
                if strand == '+':
                    start = int(a[11]) - int(a[15])
                    end = int(a[12]) + (int(a[14]) - int(a[16]))
                    position = end
                if strand == '-':
                    start = int(a[11]) - (int(a[14]) - int(a[16]))
                    end = int(a[12]) + int(a[15])
                    position = start
                adapter_dict[read_name][strand].append((adapter,
                                                        float(a[0]),
                                                        position))
    return adapter_dict

def match_index(seq, seq_to_idx):
    dist_dict, dist_list = {}, []
    # there needs to be a better/more efficient way to do this.
    for position in range(len(seq)):
        for idx_seq, idx in seq_to_idx.items():
            if idx not in dist_dict:
                dist_dict[idx] = []
            query = seq[position:position + len(idx_seq)]
            if len(query) != len(idx_seq):
                break
            else:
                dist = ld.eval(query, idx_seq)
                dist_dict[idx].append(dist)
    for idx, distances in dist_dict.items():
        dist_list.append((idx, min(distances)))
    dist_list = sorted(dist_list, key=lambda x: x[1])
    if dist_list[0][1] < 2 and dist_list[1][1] - dist_list[0][1] > 1:
        return dist_list[0][0]
    else:
        return '-'

def write_fasta_file(args, path, adapter_dict, reads, seq_to_idx, idx_to_seq):
    undirectional = args.undirectional
    barcoded = args.barcoded
    trim = args.trim

    odT = True if seq_to_idx else False

    if barcoded:
        out10X = open(path + 'R2C2_full_length_consensus_reads_10X_sequences.fasta', 'w')
    if odT:
        outdT = open(path + 'R2C2_oligodT_multiplexing.tsv', 'w')
        for idx in idx_to_seq:
            if os.path.exists(path + idx):
                shutil.rmtree(path + idx)
    else:
        out = open(path + 'R2C2_full_length_consensus_reads.fasta', 'w')
        out3 = open(path + 'R2C2_full_length_consensus_reads_left_splint.fasta', 'w')
        out5 = open(path + 'R2C2_full_length_consensus_reads_right_splint.fasta', 'w')

    for name, sequence in (tqdm(reads.items()) if args.threads==1  else reads.items()):
        adapter_plus = sorted(adapter_dict[name]['+'],
                              key=lambda x: x[2], reverse=False)
        adapter_minus = sorted(adapter_dict[name]['-'],
                              key=lambda x: x[2], reverse=False)
        plus_list_name, plus_positions = [], []
        minus_list_name, minus_positions = [], []

        for adapter in adapter_plus:
            if adapter[0] != '-':
                plus_list_name.append(adapter[0])
                plus_positions.append(adapter[2])
        for adapter in adapter_minus:
            if adapter[0] != '-':
                minus_list_name.append(adapter[0])
                minus_positions.append(adapter[2])

        if len(plus_list_name) != 1 or len(minus_list_name) != 1:
            continue
        if minus_positions[0] <= plus_positions[0]:
            continue

        if undirectional:
            direction = '+'
        elif plus_list_name[0] != minus_list_name[0]:
            if plus_list_name[0] == '5Prime_adapter':
                direction = '+'
            else:
                direction = '-'
        else:
            continue

        if odT:
            outdT.write('%s\t%s\t%s\n' %(
                name,
                mm.revcomp(sequence[minus_positions[0]-16:minus_positions[0]+4]),
                sequence[plus_positions[0]-4:plus_positions[0]+16])
            )
            reverse_index, forward_index = '-', '-'
            forward_index = match_index(sequence[plus_positions[0]-4:plus_positions[0]+16], seq_to_idx)
            reverse_index = match_index(mm.revcomp(sequence[minus_positions[0]-16:minus_positions[0]+4]), seq_to_idx)

            demux = False
            if forward_index in idx_to_seq and reverse_index not in idx_to_seq:
                direction, idx_name, demux = '-', forward_index, True
            if reverse_index in idx_to_seq and forward_index not in idx_to_seq:
                direction, idx_name, demux = '+', reverse_index, True
            if not demux:
                idx_name = 'no_index_found'

            demux_path = path + idx_name + '/'
            if not os.path.isdir(demux_path):
                os.mkdir(demux_path)

            out = open(demux_path + 'R2C2_full_length_consensus_reads.fasta', 'a+')
            out3 = open(demux_path + 'R2C2_full_length_consensus_reads_left_splint.fasta', 'a+')
            out5 = open(demux_path + 'R2C2_full_length_consensus_reads_right_splint.fasta', 'a+')

        seq = sequence[plus_positions[0]:minus_positions[0]]
        ada = sequence[max(plus_positions[0]-40, 0):minus_positions[0]+40]
        name += '_' + str(len(seq))
        if direction == '+':
            if trim:
                out.write('>%s\n%s\n' %(name, seq))
            else:
                out.write('>%s\n%s\n' %(name, ada))
            out5.write('>%s\n%s\n' %(name, mm.revcomp(sequence[:plus_positions[0]])))
            out3.write('>%s\n%s\n' %(name, sequence[minus_positions[0]:]))
            if barcoded:
                out10X.write('>%s\n%splus\n' %(name, mm.revcomp(sequence[minus_positions[0]-40:minus_positions[0]])))
        elif direction == '-':
            if trim:
                out.write('>%s\n%s\n' %(name, mm.revcomp(seq)))
            else:
                out.write('>%s\n%s\n' %(name, mm.revcomp(ada)))
            out3.write('>%s\n%s\n' %(name, mm.revcomp(sequence[:plus_positions[0]+40])))
            out5.write('>%s\n%s\n' %(name, sequence[minus_positions[0]:]))
            if barcoded:
                out10X.write('>%s\n%sminus\n' %(name, sequence[plus_positions[0]:plus_positions[0]+40]))

        if odT:
            out.close()
            out3.close()
            out5.close()

    if not odT:
        out.close()
        out3.close()
        out5.close()
    if barcoded:
        out10X.close()
    if odT:
        outdT.close()

def main(args):
    if not args.output_path.endswith('/'):
        args.output_path += '/'

    if args.config:
        progs = configReader(args.output_path, args.config)
        blat = progs['blat']
    else:
        blat = 'blat'

    if args.undirectional and args.barcoded:
        print('Error: undirectional and barcoded are mutually exclusive.')
        sys.exit(1)

    if args.threads > 1:
        num_reads = get_file_len(args.input_fasta_file)
        chunk_process(num_reads, args, blat)
    else:
        reads = read_fasta(args.input_fasta_file, False)
        
        if args.index_file:
            idx_to_seq, seq_to_idx = read_fasta(args.index_file, True)
        else:
            idx_to_seq, seq_to_idx = {}, {}

        run_blat(args.output_path, args.input_fasta_file, args.adapter_file, blat)
        adapter_dict = parse_blat(args.output_path, reads)
        write_fasta_file(args, args.output_path, adapter_dict, reads, seq_to_idx, idx_to_seq)

if __name__ == '__main__':
    args = parse_args()
    if not args.input_fasta_file or not args.adapter_file:
        print('Reads (--input_fasta_file/-i) and adapter (--adapter_file/-a) are required', file=sys.stderr)
        sys.exit(1)
    main(args)
