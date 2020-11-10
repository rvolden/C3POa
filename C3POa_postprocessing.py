#!/usr/bin/env python3
# Roger Volden and Chris Vollmers

import sys
import os
import argparse
import mappy as mm
from tqdm import tqdm
import multiprocessing as mp

def parse_args():
    '''Parses arguments.'''
    parser = argparse.ArgumentParser(description = '',
                                     add_help = True,
                                     prefix_chars = '-')
    parser.add_argument('--input_fasta_file', '-i', type=str)
    parser.add_argument('--output_path', '-o', type=str)
    parser.add_argument('--adapter_file', '-a', type=str)
    parser.add_argument('--config', '-c', type=str, action='store', default='',
                        help='If you want to use a config file to specify paths to\
                              programs, specify them here. Use for poa, racon, water,\
                              blat, and minimap2 if they are not in your path.')
    parser.add_argument('--undirectional','-u', action='store_true',
                        help='By default, your cDNA molecules are assumed to be \
                              directional with two sequences named "3Prime_adapter" \
                              and "5Prime_adapter" expected in your adapter_file in \
                              fasta format. If you add this flag your cDNA molecules \
                              are expected to be undirectional and only one sequence \
                              named "Adapter" should be in your adapter_file in fasta \
                              format')
    parser.add_argument('--trim', '-t', action='store_true',
                        help='Use this flag to trim the adapters off the ends of \
                              your sequences.')
    parser.add_argument('--barcoded', '-b', action='store_true', default=False,
                        help='Use if postprocessing 10x reads. Produces a separate \
                              file with 10x barcode sequences')
    parser.add_argument('--threads', '-n', type=int, default=1,
                        help='Number of threads to use during multiprocessing. Defaults to 1.')
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

def process(args, reads, blat, iteration):
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
    write_fasta_file(args, tmp_dir, adapter_dict, reads)

def chunk_process(num_reads, args, blat):
    '''Split the input fasta into chunks and process'''
    chunk_size = (num_reads//args.threads) + 1

    pool = mp.Pool(args.threads)
    pbar = tqdm(total=args.threads)
    iteration, current_num, tmp_reads, target = 1, 0, {}, chunk_size
    for read in mm.fastx_read(args.input_fasta_file, read_comment=False):
        tmp_reads[read[0]] = read[1]
        current_num += 1
        if current_num == target:
            pool.apply_async(process, args=(args, tmp_reads, blat, iteration), callback=lambda _: pbar.update(1))
            iteration += 1
            target = chunk_size*iteration
            if target >= num_reads:
                target = num_reads
            tmp_reads = {}
    pool.close()
    pool.join()
    pbar.close()

    os.system('cat {flc} >{flc_final}'.format(
            flc=args.output_path + 'post_tmp*/R2C2_full_length_consensus_reads.fasta',
            flc_final=args.output_path + '/R2C2_full_length_consensus_reads.fasta')
    )
    os.system('cat {flc_left} >{flc_left_final}'.format(
            flc_left=args.output_path + 'post_tmp*/R2C2_full_length_consensus_reads_left_splint.fasta',
            flc_left_final=args.output_path + '/R2C2_full_length_consensus_reads_left_splint.fasta')
    )
    os.system('cat {flc_right} >{flc_right_final}'.format(
            flc_right=args.output_path + 'post_tmp*/R2C2_full_length_consensus_reads_right_splint.fasta',
            flc_right_final=args.output_path + '/R2C2_full_length_consensus_reads_right_splint.fasta')
    )
    if args.barcoded:
        os.system('cat {flc_bc} >{flc_bc_final}'.format(
                flc_bc=args.output_path + 'post_tmp*/R2C2_full_length_consensus_reads_10X_sequences.fasta',
                flc_bc_final=args.output_path + '/R2C2_full_length_consensus_reads_10X_sequences.fasta')
        )
    os.system('rm -rf {tmps}'.format(tmps=args.output_path + 'post_tmp*'))

def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readDict = {}
    for read in mm.fastx_read(inFile, read_comment=False):
        readDict[read[0]] = read[1]
    return readDict

def run_blat(path, infile, adapter_fasta, blat):
    align_psl = path + 'adapter_to_consensus_alignment.psl'
    if not os.path.exists(align_psl) or os.stat(align_psl).st_size == 0:
        os.system('{blat} -noHead -stepSize=1 -tileSize=6 -t=DNA -q=DNA -minScore=10 \
                  -minIdentity=10 -minMatch=1 -oneOff=1 {adapters} {reads} {psl} >{blat_msgs}'\
                  .format(blat=blat, adapters=adapter_fasta, reads=infile, psl=align_psl, blat_msgs=path+'blat_msgs.log'))
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

def match_index(sequence,sequence_to_index):
    dist_dict={}
    dist_list=[]
    for position in range(0,len(sequence),1):
        for index_sequence,index in sequence_to_index.items():
            if index not in dist_dict:
                dist_dict[index]=[]
            query=sequence[position:position+len(index_sequence)]
            if len(query)!=len(index_sequence):
                break
            else:
                dist=editdistance.eval(query,index_sequence)
                dist_dict[index].append(dist)
    for index,distances in dist_dict.items():
        dist_list.append((index,min(distances)))
    dist_list=sorted(dist_list,key=lambda x: x[1])
    if dist_list[0][1]<2 and dist_list[1][1]-dist_list[0][1]>1:
        return dist_list[0][0]
    else:
        return '-'

def write_fasta_file(args, path, adapter_dict, reads):
    undirectional = args.undirectional
    barcoded = args.barcoded
    trim = args.trim

    out = open(path + 'R2C2_full_length_consensus_reads.fasta', 'w')
    out3 = open(path + 'R2C2_full_length_consensus_reads_left_splint.fasta', 'w')
    out5 = open(path + 'R2C2_full_length_consensus_reads_right_splint.fasta', 'w')
    if barcoded:
        out10X = open(path + 'R2C2_full_length_consensus_reads_10X_sequences.fasta', 'w')

    for name, sequence in (tqdm(reads.items()) if args.threads==1  else reads.items()):
        adapter_plus = sorted(adapter_dict[name]['+'],
                              key=lambda x: x[2], reverse=False)
        adapter_minus = sorted(adapter_dict[name]['-'],
                              key=lambda x: x[2], reverse=False)
        plus_list_name, plus_list_position = [], []
        minus_list_name, minus_list_position = [], []

        for adapter in adapter_plus:
            if adapter[0] != '-':
                plus_list_name.append(adapter[0])
                plus_list_position.append(adapter[2])
        for adapter in adapter_minus:
            if adapter[0] != '-':
                minus_list_name.append(adapter[0])
                minus_list_position.append(adapter[2])

        if len(plus_list_name) != 1 or len(minus_list_name) != 1:
            continue
        if minus_list_position[0] <= plus_list_position[0]:
            continue

        use = False
        if undirectional:
            direction = '+'
            use = True
        else:
            if plus_list_name[0] != minus_list_name[0]:
                use = True
                if plus_list_name[0] == '5Prime_adapter':
                    direction = '+'
                else:
                    direction = '-'
        if not use:
            continue

        seq = sequence[plus_list_position[0]:minus_list_position[0]]
        ada = sequence[max(plus_list_position[0]-40, 0):minus_list_position[0]+40]
        name += '_' + str(len(seq))
        if direction == '+':
            if trim:
                out.write('>%s\n%s\n' %(name, seq))
            else:
                out.write('>%s\n%s\n' %(name, ada))
            out5.write('>%s\n%s\n' %(name, mm.revcomp(sequence[:plus_list_position[0]])))
            out3.write('>%s\n%s\n' %(name, sequence[minus_list_position[0]:]))
            if barcoded:
                out10X.write('>%s\n%splus\n' %(name, mm.revcomp(sequence[minus_list_position[0]-40:minus_list_position[0]])))
        elif direction == '-':
            if trim:
                out.write('>%s\n%s\n' %(name, mm.revcomp(seq)))
            else:
                out.write('>%s\n%s\n' %(name, mm.revcomp(ada)))
            out3.write('>%s\n%s\n' %(name, mm.revcomp(sequence[:plus_list_position[0]+40])))
            out5.write('>%s\n%s\n' %(name, sequence[minus_list_position[0]:]))
            if barcoded:
                out10X.write('>%s\n%sminus\n' %(name, sequence[plus_list_position[0]:plus_list_position[0]+40]))
    out.close()
    out3.close()
    out5.close()
    if barcoded:
        out10X.close()

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
        reads = read_fasta(args.input_fasta_file)
        run_blat(args.output_path, args.input_fasta_file, args.adapter_file, blat)
        adapter_dict = parse_blat(args.output_path, reads)
        write_fasta_file(args, args.output_path, adapter_dict, reads)

if __name__ == '__main__':
    args = parse_args()
    main(args)
