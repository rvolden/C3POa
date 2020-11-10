#!/usr/bin/env python3
# Roger Volden

import os
import sys
import numpy as np
import argparse
import multiprocessing as mp
import mappy as mm
from conk import conk
from tqdm import tqdm

path = '/'.join(os.path.realpath(__file__).split('/')[:-1])+'/bin/'
sys.path.append(os.path.abspath(path))

from preprocess import preprocess
from call_peaks import call_peaks
from determine_consensus import determine_consensus

def parse_args():
    '''Parses arguments.'''
    parser = argparse.ArgumentParser(description = 'Makes consensus sequences \
                                                    from R2C2 reads.',
                                     add_help = True,
                                     prefix_chars = '-')
    required = parser.add_argument_group('required arguments')
    required.add_argument('--reads', '-r', type=str, action='store', required=True,
                          help='FASTQ file that contains the long R2C2 reads.')
    required.add_argument('--splint_file', '-s', type=str, action='store', required=True,
                        help='Path to the splint FASTA file.')
    parser.add_argument('--out_path', '-o', type=str, action='store', default=os.getcwd(),
                        help='''Directory where all the files will end up.
                                Defaults to your current directory.''')
    parser.add_argument('--config', '-c', type=str, action='store', default='',
                        help='''If you want to use a config file to specify paths to
                                programs, specify them here. Use for racon and blat
                                if they are not in your path.''')
    parser.add_argument('--lencutoff', '-l', type=int, action='store', default=1000,
                        help='''Sets the length cutoff for your raw sequences. Anything
                                shorter than the cutoff will be excluded. Defaults to 1000.''')
    parser.add_argument('--mdistcutoff', '-d', type=int, action='store', default=500,
                        help='''Sets the median distance cutoff for consensus sequences.
                                Anything shorter will be excluded. Defaults to 500.''')
    parser.add_argument('--zero', '-z', action='store_false', default=True,
                        help='Use to exclude zero repeat reads. Defaults to True (includes zero repeats).')
    parser.add_argument('--numThreads', '-n', type=int, default=1,
                        help='Number of threads to use during multiprocessing. Defaults to 1.')
    parser.add_argument('--groupSize', '-g', type=int, default=1000,
                        help='Number of reads processed by each thread in each iteration. Defaults to 1000.')
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

def rounding(x, base):
    '''Rounds to the nearest base, we use 50'''
    return int(base * round(float(x)/base))

def analyze_reads(args, reads, splint_dict, adapter_dict, adapter_set, iteration, racon):
    penalty, iters, window, order = 20, 3, 41, 2
    for read in reads:
        name, seq, qual = read[0], read[1], read[2]
        seq_len = len(seq)
        if not adapter_dict.get(name):
            continue
        strand = adapter_dict[name][1]
        if strand == '-':
            # use reverse complement of the splint
            splint = splint_dict[adapter_dict[name][0]][1]
        else:
            splint = splint_dict[adapter_dict[name][0]][0]
        scores = conk.conk(splint, seq, penalty)
        peaks = call_peaks(scores, args.mdistcutoff, iters, window, order)
        if not list(peaks):
            continue
        peaks = list(peaks + len(splint)//2)
        for i in range(len(peaks)-1, -1, -1):
            if peaks[i] >= seq_len:
                del peaks[i]
        if not peaks:
            continue

        # check for outliers in subread length
        subreads, qual_subreads = [], []
        if len(peaks) > 1:
            subread_lens = np.diff(peaks)
            subread_lens = [rounding(x, 50) for x in subread_lens]
            median_subread_len = np.median(subread_lens)
            for i in range(len(subread_lens)):
                bounds = [peaks[i], peaks[i+1]]
                if median_subread_len*0.8 <= subread_lens[i] <= median_subread_len*1.2:
                    subreads.append(seq[bounds[0]:bounds[1]])
                    qual_subreads.append(qual[bounds[0]:bounds[1]])
            if peaks[0] > 100:
                subreads.append(seq[:peaks[0]])
                qual_subreads.append(qual[:peaks[0]])
            if seq_len - peaks[-1] > 100:
                subreads.append(seq[peaks[-1]:])
                qual_subreads.append(qual[peaks[-1]:])
        else:
            subreads.append(seq[:peaks[0]])
            qual_subreads.append(qual[:peaks[0]])
            subreads.append(seq[peaks[0]:])
            qual_subreads.append(qual[peaks[0]:])

        tmp_dir = args.out_path + adapter_dict[name][0] + '/tmp' + str(iteration) + '/'
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
        subread_file = tmp_dir + 'subreads.fastq'

        consensus, repeats = determine_consensus(
            args, read, subreads, qual_subreads,
            racon, tmp_dir, subread_file
        )

        if consensus:
            avg_qual = round(sum([ord(x)-33 for x in qual])/seq_len, 2)
            cons_len = len(consensus)
            final_out = open(tmp_dir + '/R2C2_Consensus.fasta', 'a+')
            print('>' + name + '_' + '_'.join([str(x) for x in [avg_qual, seq_len, repeats, cons_len]]), file=final_out)
            print(consensus, file=final_out)
            final_out.close()

def main(args):
    if not args.out_path.endswith('/'):
        args.out_path += '/'
    if not os.path.exists(args.out_path):
        os.mkdir(args.out_path)
    log_file = open(args.out_path + 'c3poa.log', 'w+')

    if args.config:
        progs = configReader(args.out_path, args.config)
        racon = progs['racon']
        blat = progs['blat']
    else:
        racon = 'racon'
        blat = 'blat'

    tmp_dir = args.out_path + 'tmp/'
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)

    # read in the file and preprocess
    read_list, total_reads = [], 0
    short_reads = 0
    tmp_fasta = tmp_dir + 'R2C2_temp_for_BLAT.fasta'
    align_psl = tmp_dir + 'splint_to_read_alignments.psl'

    '''
    Workflow for blat parallelization:
        If the psl file exists, don't parallelize and set make_tmp_fasta to false

        else:
            use the same chunking method from postprocessing to make the input
            fasta files for blat (just like what we're doing now), but the chunks
            will be bigger.

            Implement the mp pool in preprocess.py
            after, cat all of the output psl files into a master psl file
            from there, read in that master psl (what's already implemented)
    '''

    make_tmp_fasta, tmp_adapter_dict = False, {}
    if not os.path.exists(align_psl) or os.stat(align_psl).st_size == 0:
        make_tmp_fasta = True
        tmp_fa_fh = open(tmp_fasta, 'w')
    for read in mm.fastx_read(args.reads, read_comment=False):
        if len(read[1]) < args.lencutoff:
            short_reads += 1
            continue
        read_list.append(read) # [(read_id, seq, qual), ...]
        if make_tmp_fasta:
            print('>' + read[0] + '\n' + read[1], file=tmp_fa_fh)
        tmp_adapter_dict[read[0]] = [[None, 1, None]] # [adapter, matches, strand]
        total_reads += 1
    if make_tmp_fasta:
        tmp_fa_fh.close()

    adapter_dict, adapter_set, no_splint = preprocess(blat, args.out_path, tmp_dir, read_list, args.splint_file, tmp_adapter_dict)
    for adapter in adapter_set:
        if not os.path.exists(args.out_path + adapter):
            os.mkdir(args.out_path + adapter)

    print('No splint reads:', no_splint, file=log_file)
    print('Under len cutoff:', short_reads, file=log_file)
    print('Total thrown away reads:', short_reads + no_splint, file=log_file)
    log_file.close()
 
    splint_dict = {}
    for splint in mm.fastx_read(args.splint_file, read_comment=False):
        splint_dict[splint[0]] = [splint[1]]
        splint_dict[splint[0]].append(mm.revcomp(splint[1]))

    pool = mp.Pool(args.numThreads)
    iteration = 1
    pbar = tqdm(total=total_reads//args.groupSize+1)
    for step in range(0, total_reads, args.groupSize):
        interval1 = step
        interval2 = min(total_reads, step+args.groupSize)

        pool.apply_async(analyze_reads, \
            args=(args, read_list[interval1:interval2], splint_dict, adapter_dict, adapter_set, iteration, racon,), \
            callback=lambda _: pbar.update(1) \
        )
        iteration += 1
    pool.close()
    pool.join()
    pbar.close()

    for adapter in adapter_set:
        os.system('cat {adapter_dir} >{cons}'.format(
                adapter_dir=args.out_path + adapter + '/tmp*/R2C2_Consensus.fasta',
                cons=args.out_path + adapter + '/R2C2_Consensus.fasta')
        )
        os.system('cat {adapter_dir} >{subreads}'.format(
                adapter_dir=args.out_path + adapter + '/tmp*/subreads.fastq',
                subreads=args.out_path + adapter + '/R2C2_Subreads.fastq')
        )
        os.system('rm -rf {tmps}'.format(tmps=args.out_path + adapter + '/tmp*'))

if __name__ == '__main__':
    args = parse_args()
    mp.set_start_method("spawn")
    main(args)
