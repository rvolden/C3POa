#!/usr/bin/env python3
# Roger Volden

import pyabpoa as poa
import mappy as mm
import os
import subprocess
from consensus import pairwise_consensus

def determine_consensus(args, read, subreads, sub_qual, racon, tmp_dir, subread_file):
    name, seq, qual = read[0], read[1], read[2]
    repeats = len(subreads)

    if repeats == 2 and args.zero:
        final_cons = zero_repeats(name, seq, qual, subreads, sub_qual, subread_file)
        if final_cons and len(final_cons) >= args.mdistcutoff:
            return final_cons, 0

    # subread is the master subread fastq for this group
    subread_fh = open(subread_file, 'a+')
    # overlap file is where the mappy alignment will go (req. by racon)
    overlap_file = tmp_dir + '{name}_overlaps.paf'.format(name=name)
    overlap_fh = open(overlap_file, 'w+')
    # temporary subreads specific for the current read (req. by racon)
    tmp_subread_file = tmp_dir + '{name}_subreads.fastq'.format(name=name)
    tmp_subread_fh = open(tmp_subread_file, 'w+')

    # align subreads together using abPOA
    poa_aligner = poa.msa_aligner(match=5)
    if repeats == 2:
        res = poa_aligner.msa(subreads, out_cons=False, out_msa=True)
        if not res.msa_seq:
            subread_fh.close()
            overlap_fh.close()
            tmp_subread_file.close()
            os.system('rm {tmp_files}'.format(tmp_files=' '.join([tmp_subread_file, overlap_file])))
            return '', 0
        abpoa_cons = pairwise_consensus(res.msa_seq, subreads, sub_qual)
    else:
        res = poa_aligner.msa(subreads, out_cons=True, out_msa=True)
        if not res.cons_seq:
            os.system('rm {tmp_files}'.format(tmp_files=' '.join([tmp_subread_file, overlap_file])))
            return '', 0
        abpoa_cons = res.cons_seq[0]

    # have to write out the consensus seq because it's going to get polished by racon
    abpoa_fasta = tmp_dir + '{name}_abpoa.fasta'.format(name=name)
    abpoa_fasta_fh = open(abpoa_fasta, 'w+')
    print('>{name}\n{seq}\n'.format(name=name, seq=abpoa_cons), file=abpoa_fasta_fh)
    abpoa_fasta_fh.close()

    # map each of the subreads to the poa consensus
    mm_align = mm.Aligner(seq=abpoa_cons, preset='map-ont')
    for i in range(repeats):
        subread = subreads[i]
        q = sub_qual[i]
        for hit in mm_align.map(subread):
            qname = name + '_subread_' + str(i)
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                qname, str(len(subread)), hit.q_st, hit.q_en,
                hit.strand, name, hit.ctg_len, hit.r_st,
                hit.r_en, hit.mlen, hit.blen, hit.mapq), file=overlap_fh)
            print('@{name}\n{sub}\n+\n{q}'.format(name=qname, sub=subread, q=q), file=subread_fh)
            print('@{name}\n{sub}\n+\n{q}'.format(name=qname, sub=subread, q=q), file=tmp_subread_fh)
    subread_fh.close()
    overlap_fh.close()
    tmp_subread_fh.close()

    racon_cons_file = tmp_dir + '{name}_racon_cons.fasta'.format(name=name)
    racon_cons_fh = open(racon_cons_file, 'w+')
    racon_msgs_fh = open(tmp_dir + 'racon_messages.log', 'w+')

    # polish poa cons with the subreads
    subprocess.run([racon, tmp_subread_file, overlap_file, abpoa_fasta, '-q', '5', '-t', '1'],
                   stdout=racon_cons_fh, stderr=racon_msgs_fh)
    racon_cons_fh.close()
    racon_msgs_fh.close()

    final_cons = ''
    for read in mm.fastx_read(racon_cons_file, read_comment=False):
        final_cons = read[1]

    tmp_files = ' '.join([overlap_file, tmp_subread_file, abpoa_fasta, racon_cons_file])
    os.system('rm {tmp_files}'.format(tmp_files=tmp_files))

    return final_cons, repeats

def zero_repeats(name, seq, qual, subreads, sub_qual, subread_file):
    # subread is the master subread fastq for this group
    subread_fh = open(subread_file, 'a+')
    for i in range(len(subreads)):
        print('@{name}\n{sub}\n+\n{q}'.format(name=name + '_subread_' + str(i),
                                              sub=subreads[i],
                                              q=sub_qual[i]),
                                              file=subread_fh)
    subread_fh.close()

    mappy_res = []
    mm_align = mm.Aligner(seq=subreads[0], preset='map-ont', scoring=(20, 7, 10, 5))
    for hit in mm_align.map(subreads[1]):
        mappy_res = [hit.r_st, hit.r_en, hit.q_st, hit.q_en]
    if not mappy_res:
        return ''

    left = subreads[0][:mappy_res[0]]
    right = subreads[1][mappy_res[3]:]
    overlap_seq1 = subreads[0][mappy_res[0]:mappy_res[1]]
    overlap_qual1 = sub_qual[0][mappy_res[0]:mappy_res[1]]
    overlap_seq2 = subreads[1][mappy_res[2]:mappy_res[3]]
    overlap_qual2 = sub_qual[1][mappy_res[2]:mappy_res[3]]

    poa_aligner = poa.msa_aligner(match=5)
    res = poa_aligner.msa([overlap_seq1, overlap_seq2], out_cons=False, out_msa=True)
    if not res.msa_seq:
        return ''
    abpoa_cons = pairwise_consensus(res.msa_seq, [overlap_seq1, overlap_seq2], [overlap_qual1, overlap_qual2])
    corrected_cons = left + abpoa_cons + right
    return corrected_cons
