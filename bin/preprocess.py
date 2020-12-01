#!/usr/bin/env python3
# Roger Volden

import os
import sys
import mappy as mm
from tqdm import tqdm
import multiprocessing as mp

def preprocess(blat, args, tmp_dir, tmp_adapter_dict, num_reads):
    tmp_fasta = tmp_dir + 'R2C2_temp_for_BLAT.fasta'
    align_psl = tmp_dir + 'splint_to_read_alignments.psl'

    # skip the alignment if the psl file already exists
    if not os.path.exists(align_psl) or os.stat(align_psl).st_size == 0:
        print('Aligning splints to reads with blat', file=sys.stderr)
        chunk_process(num_reads, args, blat)
    else:
        print('Reading existing psl file', file=sys.stderr)

    adapter_set = set()
    with open(align_psl) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            line = line.split('\t')
            read_name, adapter, strand = line[9], line[13], line[8]
            gaps, score = float(line[5]), float(line[0])
            if gaps < 50 and score > 50:
                tmp_adapter_dict[read_name].append([adapter, float(line[0]), strand])
                adapter_set.add(adapter)

    adapter_dict = {} # read_id: [adapter, strand]
    no_splint_reads = 0
    for name, alignments in tmp_adapter_dict.items():
        best = sorted(alignments, key=lambda x: x[1], reverse=True)[0]
        if not best[0]:
            no_splint_reads += 1
            continue
        adapter_set.add(best[0])
        adapter_dict[name] = [best[0], best[2]]
    return adapter_dict, adapter_set, no_splint_reads

def process(args, reads, blat, iteration):
    tmp_dir = args.out_path + 'pre_tmp_' + str(iteration) + '/'
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)
    tmp_fa = tmp_dir + 'tmp_for_blat.fasta'
    tmp_fa_fh = open(tmp_fa, 'w+')
    for header, seq in reads.items():
        print('>' + header, file=tmp_fa_fh)
        print(seq, file=tmp_fa_fh)
    tmp_fa_fh.close()
    align_psl = tmp_dir + 'tmp_splint_aln.psl'
    b_msgs = tmp_dir + 'blat_messages.log'

    os.system('{blat} -noHead -stepSize=1 -t=DNA -q=DNA -minScore=15 \
              -minIdentity=10 {splint} {reads} {psl} >{blat_msgs}'
              .format(blat=blat, splint=args.splint_file, reads=tmp_fa, psl=align_psl, blat_msgs=b_msgs))
    os.remove(tmp_fa)

def chunk_process(num_reads, args, blat):
    '''Split the input fasta into chunks and process'''
    if args.blatThreads:
        chunk_size = (num_reads // args.numThreads) + 1
    else:
        chunk_size = args.groupSize
    if chunk_size > num_reads:
        chunk_size = num_reads

    pool = mp.Pool(args.numThreads)
    pbar = tqdm(total=num_reads // chunk_size + 1, desc='Preprocessing')
    iteration, current_num, tmp_reads, target = 1, 0, {}, chunk_size
    for read in mm.fastx_read(args.reads, read_comment=False):
        if len(read[1]) < args.lencutoff:
            continue
        tmp_reads[read[0]] = read[1]
        current_num += 1
        if current_num == target:
            pool.apply_async(process, args=(args, tmp_reads, blat, iteration), callback=lambda _: pbar.update(1))
            iteration += 1
            target = chunk_size * iteration
            if target >= num_reads:
                target = num_reads
            tmp_reads = {}
    pool.close()
    pool.join()
    pbar.close()

    os.system('cat {psls} >{psl_final}'.format(
            psls=args.out_path + 'pre_tmp_*/tmp_splint_aln.psl',
            psl_final=args.out_path + 'tmp/splint_to_read_alignments.psl')
    )
    os.system('rm -rf {tmps}'.format(tmps=args.out_path + 'pre_tmp*'))
