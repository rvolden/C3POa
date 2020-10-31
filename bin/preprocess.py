#!/usr/bin/env python3
# Roger Volden

import os
import sys

def preprocess(blat, out_dir, tmp_dir, reads, splint_file, tmp_adapter_dict):
    tmp_fasta = tmp_dir + 'R2C2_temp_for_BLAT.fasta'
    align_psl = tmp_dir + 'splint_to_read_alignments.psl'

    # skip the alignment if the psl file already exists
    if not os.path.exists(align_psl) or os.stat(align_psl).st_size == 0:
        print('Aligning splints to reads with blat', file=sys.stderr)

        # run blat on all reads to make the psl file
        os.system('{blat} -noHead -stepSize=1 -t=DNA -q=DNA -minScore=15 \
                  -minIdentity=10 {splint} {reads} {psl}'\
                  .format(blat=blat, splint=splint_file, reads=tmp_fasta, psl=align_psl))
        os.remove(tmp_fasta)
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
    for read in reads:
        name, seq, qual = read[0], read[1], read[2]
        # get the alignments with the most matches
        best = sorted(tmp_adapter_dict[name], key=lambda x:x[1], reverse=True)[0]
        if not best[0]:
            no_splint_reads += 1
            continue
        adapter_set.add(best[0])
        adapter_dict[read[0]] = [best[0], best[2]]
    return adapter_dict, adapter_set, no_splint_reads
