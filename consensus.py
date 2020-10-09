#!/usr/bin/env python3
# Roger Volden

'''
Take a FASTA file with aligned reads and a FASTQ file to make a consensus
based on quality as well as base frequency at certain positions.
FASTA file with aligned reads can also come as .fasta from EMBOSS needle.

Usage: python3 consensus.py aligned.fastq reads.fastq >consensus.fasta
'''

import sys

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
    print('>consensus')
    print(consensus.replace('-', ''))

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

def fastaReader(inFile):
    '''Reads in FASTA files.'''
    tempSeqs, headers, sequences = [], [], []
    for line in inFile:
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            headers.append(line.split()[0][1:])
        # covers the case where the file ends while reading sequences
        if line.startswith('>'):
            sequences.append(''.join(tempSeqs).upper())
            tempSeqs = []
        else:
            tempSeqs.append(line)
    sequences.append(''.join(tempSeqs).upper())
    return sequences[1:]

def fastqReader(inFile):
    '''Reads in FASTQ files. Only returns sequence and quality lines.'''
    sequences, quality, lineNum = [], [], 0
#    print('fastq')
    for line in inFile:
        line = line.rstrip()
#        print(line)
        if not line:
            continue
        # sequences
        if lineNum % 4 == 1:
            sequences.append(line)
        # quality lines
        elif lineNum % 4 == 3:
            quality.append(line)
        lineNum += 1
    return sequences, quality

def main():
    '''
    Reads in files and organizes information.
    Calls consensus, which actually builds the consensus sequence.
    '''
    alignedSeqs = fastaReader(open(sys.argv[1]))
    fastqSeqs, qualityScores = fastqReader(open(sys.argv[2]))
    seqDict = {} # seq:quality
    for i in range(len(fastqSeqs)):
        seqDict[fastqSeqs[i]] = qualityScores[i]
    consensus(alignedSeqs, seqDict)

main()
