#!/usr/bin/env python3
# Roger Volden

'''
Take a FASTA file with aligned reads and a FASTQ file to make a consensus
based on quality as well as base frequency at certain positions.
FASTA file with aligned reads can also come as .fasta from EMBOSS needle.

Usage: python3 consensus.py aligned.fastq reads.fastq >consensus.fasta
'''

def consensus(sequences, qualityDict):
    '''
    Makes a consensus sequence based on base frequency and quality.
    sequences: list of aligned sequences.
    qualityDict: dictionary of sequences : quality scores.
    Returns a consensus sequence.
    '''
    consensus = ''
    seqA, seqB = sequences[0], sequences[1]
    seqAq, seqBq = qualityDict[seqA.replace('-', '')], qualityDict[seqB.replace('-', '')]
    seqAqual = normalizeLen(seqA, seqAq)
    seqBqual = normalizeLen(seqB, seqBq)

    i = 0
    while i != len(seqA): # iterate by position
        if seqA[i] == seqB[i]: # match
            consensus += seqA[i]
        if seqA[i] != seqB[i] and seqA[i] != '-' and seqB[i] != '-': # mismatch
            if ord(seqAqual[i]) > ord(seqBqual[i]):
                consensus += seqA[i]
            else:
                consensus += seqB[i]
        if seqA[i] == '-' or seqB[i] == '-': # gap to bases
            gapLen = 1 # where to start gap chunk
            if seqA[i] == '-': # which seq to check
                gapSeq = seqA
            else:
                gapSeq = seqB
            try:
                while gapSeq[i + gapLen] == '-': # extend to length of gap
                    gapLen += 1
            except IndexError: # if gap at end
                gapLen = 1
            if avgQual(seqAqual, i, gapLen) > avgQual(seqBqual, i, gapLen):
                consensus += seqA[i:i+gapLen]
            else:
                consensus += seqB[i:i+gapLen]
            i += gapLen
            continue
        i += 1
    return consensus.replace('-', '')

def avgQual(qual, i, gapLen):
    '''Returns average quality of a segment.'''
    return sum(ord(x) for x in list(qual[i:i+gapLen]))/gapLen

def normalizeLen(seq, quality):
    '''
    Inserts avg quality scores based on surrounding quality scores
    where there are gaps in the sequence.
    Returns a new quality string that's the same len as the sequence.
    '''
    seqIndex, qualIndex = 0, 0
    newQuality = ''
    while qualIndex < len(quality):
        if seq[seqIndex] != '-':
            newQuality += quality[qualIndex]
            qualIndex += 1
            seqIndex += 1
        elif seq[seqIndex] == '-' and qualIndex == 0:
            newQuality += quality[qualIndex]
            seqIndex += 1
        else:
            newQuality += chr(int((ord(quality[qualIndex-1]) + ord(quality[qualIndex]))/2))
            seqIndex += 1
    if len(seq) != len(newQuality):
        gapLen = 0
        while seq[-1-gapLen] == '-':
            newQuality += newQuality[-1]
            gapLen += 1
    return newQuality

def pairwise_consensus(poa_subreads, subreads, sub_quals):
    seqDict = {}
    for i in range(len(subreads)):
        seqDict[subreads[i]] = sub_quals[i]
    pw_cons = consensus(poa_subreads, seqDict)
    return pw_cons

