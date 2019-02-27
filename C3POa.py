#!/usr/bin/env python3
# Roger Volden and Chris Vollmers
# Last updated: 27 Feb 2019

'''
Concatemeric Consensus Caller with Partial Order Alignments (C3POa)

Analyses reads by reading them in, doing self-self alignments, calling
peaks in alignment scores,  splitting reads, aligning those to each other,
and giving back a consensus sequence.

Usage:
    python3 C3POa.py --reads reads.fastq [--path /current/directory]

Dependencies:
    Python 3.6
    NumPy 1.13.3
    poa v1.0.0 Revision: 1.2.2.9500
    gonk
    minimap2 2.7-r654
    racon

02/08/2018 Release note:
    By default, this will now output what we call zero repeat reads along
    with the rest of the R2C2 reads. Zero repeat reads are reads that contain
    a splint with incomplete portions of your original molecule on each side.
    If there's an overlap, it'll align the portions that overlap and
    concatenate the rest of the read together to try and make a contiguous
    read. These reads are very similar to normal 1D reads, but there are a few
    cases where there is a slight improvement. There will be an option to
    remove these reads in postprocessing.

02/27/2019 Release note:
    I have changed the aligner to gonk over water for faster alignments. Because
    gonk does not do a complete alignment, I have removed support for zero repeat
    reads. This is because zero repeat reads only increase the number of reads
    you end up with instead of increasing the quality of the dataset. The old
    version of C3POa can be found at https://github.com/rvolden/C3POa/tree/water.
'''

import os
import sys
import numpy as np
import argparse
from time import time

def argParser():
    '''Parses arguments.'''
    parser = argparse.ArgumentParser(description = 'Makes consensus sequences \
                                                    from R2C2 reads.',
                                     add_help = True,
                                     prefix_chars = '-')
    required = parser.add_argument_group('required arguments')
    required.add_argument('--reads', '-r', type=str, action='store', required=True,
                          help='FASTQ file that contains the long R2C2 reads.')
    parser.add_argument('--path', '-p', type=str, action='store', default=os.getcwd(),
                        help='Directory where all the files are/where they will end up.\
                              Defaults to your current directory.')
    parser.add_argument('--matrix', '-m', type=str, action='store',
                        default='NUC.4.4.mat',
                        help='Score matrix to use for poa.\
                              Defaults to NUC.4.4.mat.')
    parser.add_argument('--config', '-c', type=str, action='store', default='',
                        help='If you want to use a config file to specify paths to\
                              programs, specify them here. Use for poa, racon, gonk,\
                              blat, and minimap2 if they are not in your path.')
    parser.add_argument('--slencutoff', '-l', type=int, action='store', default=1000,
                        help='Sets the length cutoff for your raw sequences. Anything\
                              shorter than the cutoff will be excluded. Defaults to 1000.')
    parser.add_argument('--mdistcutoff', '-d', type=int, action='store', default=500,
                        help='Sets the median distance cutoff for consensus sequences.\
                              Anything shorter will be excluded. Defaults to 500.')
    parser.add_argument('--output', '-o', type=str, action='store',
                        default='R2C2_Consensus.fasta',
                        help='FASTA file that the consensus gets written to.\
                              Defaults to R2C2_Consensus.fasta.')
    parser.add_argument('--timer', '-t', action='store_true', default=False,
                        help='Prints how long each dependency takes to run.\
                              Defaults to False.')
    parser.add_argument('--figure', '-f', action='store_true', default=False,
                        help='Use if you want to output a histogram of scores.')
    return vars(parser.parse_args())

def configReader(configIn):
    '''Parses the config file.'''
    progs = {}
    for line in open(configIn):
        if line.startswith('#') or not line.rstrip().split():
            continue
        line = line.rstrip().split('\t')
        progs[line[0]] = line[1]
    # should have minimap, poa, racon, gonk, consensus
    # check for extra programs that shouldn't be there
    possible = set(['poa', 'minimap2', 'gonk', 'consensus', 'racon', 'blat'])
    inConfig = set()
    for key in progs.keys():
        inConfig.add(key)
        if key not in possible:
            raise Exception('Check config file')
    # check for missing programs
    # if missing, default to path
    for missing in possible-inConfig:
        if missing == 'consensus':
            path = 'consensus.py'
        else:
            path = missing
        progs[missing] = path
        sys.stderr.write('Using ' + str(missing)
                         + ' from your path, not the config file.\n')
    return progs

args = argParser()
if args['config']:
    progs = configReader(args['config'])
    minimap2 = progs['minimap2']
    poa = progs['poa']
    racon = progs['racon']
    gonk = progs['gonk']
    consensus = progs['consensus']
else:
    minimap2, poa, racon, gonk = 'minimap2', 'poa', 'racon', 'gonk'
    consensus = 'consensus.py'

consensus = 'python3 ' + consensus
path = args['path'] + '/'
temp_folder = path + '/' + 'tmp1'
input_file = args['reads']
score_matrix = args['matrix']

seqLenCutoff = args['slencutoff']
medDistCutoff = args['mdistcutoff']

out_file = args['output']
timer = args['timer']
figure = args['figure']
subread_file = 'subreads.fastq'
os.chdir(path)
sub = open(path + '/' + subread_file, 'w')
os.system('rm -r ' + temp_folder)
os.system('mkdir ' + temp_folder)

def revComp(sequence):
    '''Returns the reverse complement of a sequence'''
    bases = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', '-':'-'}
    return ''.join([bases[x] for x in list(sequence)])[::-1]

def split_read(split_list, sequence, out_file1, qual, out_file1q, name):
    '''
    split_list : list, peak positions
    sequence : str
    out_file1 : output FASTA file
    qual : str, quality line from FASTQ
    out_file1q : output FASTQ file
    name : str, read ID

    Writes sequences to FASTA and FASTQ files.
    Returns number of repeats in the sequence.
    '''
    out_F = open(out_file1, 'w')
    out_Fq = open(out_file1q, 'w')
    for i in range(len(split_list) - 1):
        split1 = split_list[i]
        split2 = split_list[i+1]
        if len(sequence[split1:split2]) > 30:
            out_F.write('>' + str(i + 1) + '\n' \
                        + sequence[split1:split2] + '\n')
            out_Fq.write('@' + str(i + 1) + '\n' \
                         + sequence[split1:split2] + '\n+\n' \
                         + qual[split1:split2] + '\n')
            sub.write('@' + name + '_' + str(i + 1) +' \n' \
                      + sequence[split1:split2] + '\n+\n' \
                      + qual[split1:split2] + '\n')

    if len(sequence[:split_list[0]]) > 50:
        out_Fq.write('@' + str(0) + '\n' \
                     + sequence[0:split_list[0]] + '\n+\n' \
                     + qual[0:split_list[0]] + '\n')
        sub.write('@' + name + '_' + str(0) + '\n' \
                  + sequence[0:split_list[0]] + '\n+\n' \
                  + qual[0:split_list[0]] + '\n')

    if len(sequence[split2:]) > 50:
        out_Fq.write('@' + str(i + 2) + '\n' \
                     + sequence[split2:] + '\n+\n' \
                     + qual[split2:] + '\n')
        sub.write('@' + name + '_' + str(i + 2) + '\n' \
                  + sequence[split2:] + '\n+\n' \
                  + qual[split2:] + '\n')
    repeats = str(int(i + 1))
    out_F.close()
    out_Fq.close()
    return repeats

def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readDict = {}
    for line in open(inFile):
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            readDict[line[1:]] = ''
            lastHead = line[1:]
        else:
            readDict[lastHead] += line
    return readDict

def rounding(x, base):
    '''Rounds to the nearest base, we use 50'''
    return int(base * round(float(x)/base))

def makeFig(scoreList, peaks, seed, filtered_peaks):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mplpatches
    plt.style.use('BME163')
    plt.figure(figsize = (10, 5))
    hist = plt.axes([0.1, 0.1, 8/10, 4/5], frameon = True)

    xlist = [x for x in range(0, len(filtered_peaks))]
    hist.plot(xlist, filtered_peaks, color =  (0, 68/255, 85/255),
              lw = 1, zorder = 550)
    ylim = max(scoreList) * 1.1
    ymin = min(filtered_peaks)*1.5
    xlim = len(scoreList)

    for i in range(len(scoreList)):
        if np.in1d(i, peaks):
            color = (0.96, 0.43, 0.2)
            peakMark = mplpatches.Rectangle((i-12.5, filtered_peaks[i]), 25, ylim,
                                            lw=0, fc=color, zorder=0)
            hist.add_patch(peakMark)
        color = (0, 191/255, 165/255)
        score = mplpatches.Rectangle((i, 0), 1, scoreList[i],
                                     lw=0, fc=color, zorder=100)
        hist.add_patch(score)

    hist.set_ylim(ymin, ylim)
    hist.set_xlim(0, xlim)
    hist.set_ylabel('Alignment Score', fontsize = 11, labelpad = 6.5)
    hist.set_xlabel('Read position', fontsize = 11, labelpad = 6)
    hist.tick_params(axis='both',which='both',\
                     bottom='on', labelbottom='on',\
                     left='on', labelleft='on',\
                     right='off', labelright='off',\
                     top='off', labeltop='off')

    plt.savefig('plumetest.png', dpi = 600)
    plt.close()
    sys.exit()

def savitzky_golay(y, window_size, order, deriv=0, rate=1, returnScoreList=False):
    '''
    Smooths over data using a Savitzky Golay filter
    This can either return a list of scores, or a list of peaks

    y : array-like, score list
    window_size : int, how big of a window to smooth
    order : what order polynomial
    returnScoreList : bool
    '''
    from math import factorial
    y = np.array(y)
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order + 1)
    half = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half, half + 1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    filtered = np.convolve( m[::-1], y, mode='valid')

    if returnScoreList:
        return np.convolve( m[::-1], y, mode='valid')

    # set everything between 1 and -inf to 1
    posFiltered = []
    for i in range(len(filtered)):
        if 1 > filtered[i] >= -np.inf:
            posFiltered.append(1)
        else:
            posFiltered.append(filtered[i])

    # use slopes to determine peaks
    peaks = []
    slopes = np.diff(posFiltered)
    la = 45 # how far in sequence to look ahead
    for i in range(len(slopes) - 50):
        if i > len(slopes) - la: # probably irrelevant now
            dec = all(slopes[i+x] <= 0 for x in range(1, 50))
            if slopes[i] > 0 and dec:
                if i not in peaks:
                    peaks.append(i)
        else:
            dec = all(slopes[i+x] <= 0 for x in range(1, la))
            if slopes[i] > 0 and dec:
                peaks.append(i)
    return peaks

def callPeaks(scoreList):
    '''
    scoreList : list of scores
    returns a sorted list of all peaks
    '''
    maxScore = max(scoreList)
    noise = maxScore*0.05

    for j in range(len(scoreList)):
        if scoreList[j] <= noise:
            scoreList[j] = 1

    # Smooth over the data
    smoothedScores = savitzky_golay(scoreList, 21, 2, deriv = 0,
                                     rate = 1, returnScoreList = True)
    peaks = savitzky_golay(smoothedScores, 51, 1, deriv = 0,
                            rate = 1, returnScoreList = False)
    # Add all of the smoothed peaks to list of all peaks
    sortedPeaks = sorted(list(set(peaks)))

    finalPeaks = [sortedPeaks[0]]
    for i in range(1, len(sortedPeaks)):
        if sortedPeaks[i-1] < sortedPeaks[i] < sortedPeaks[i-1] + 100:
            continue
        else:
            finalPeaks.append(sortedPeaks[i])
    if figure:
        return finalPeaks, smoothedScores

    # calculates the median distance between detected peaks
    forMedian = []
    for i in range(len(finalPeaks) - 1):
        forMedian.append(finalPeaks[i+1] - finalPeaks[i])
    forMedian = [rounding(x, 50) for x in forMedian]
    medianDistance = np.median(forMedian)
    return finalPeaks, medianDistance

def parse_file(scores):
    '''
    scores : gonk output file
    Returns:
        scoreList : list, diagonal alignment scores
    '''
    scoreList = []
    for line in open(scores):
        line = line.rstrip().split(':')
        value = int(line[1])
        scoreList.append(value)
    return scoreList

def runGonk(seq1, seq2):
    '''Runs gonk using the sequences given by split_SW'''
    go_start = time()
    os.system('{0} -a seq1.fasta -b seq2.fasta -p 20 &>>gonk_messages'.format(gonk))
    go_stop = time()
    if timer:
        print('gonk took ' + str(go_stop - go_start) + ' seconds to run.')
    scores = 'SW_PARSE.txt'
    scoreList = parse_file(scores)
    os.system('rm {0}'.format(scores))
    return scoreList

def split_SW(name, seed, seq):
    '''
    Takes a sequence and does the gonk alignment to itself
    Returns a list of scores from summing diagonals from the
    alignment matrix.
    name (str): the sequence header
    seq (str): nucleotide sequence
    '''
    total = len(seq)
    reverse = False
    if seed + 1000 > total:
        start = max(0, seed-1000)
        seq1 = revComp(seq[start:seed])
        seq = revComp(seq)
        reverse = True
    else:
        seq1 = seq[seed:seed+1000]

    align_file1 = open('seq1.fasta', 'w')
    align_file1.write('>' + name + '\n' + seq1 + '\n')
    align_file1.close()
    align_file2 = open('seq2.fasta', 'w')
    align_file2.write('>' + name + '\n')
    for i in range(0, len(seq), 5000):
        align_file2.write(seq[i:i+5000] + '\n')
    align_file2.close()

    scoreList = runGonk(seq1, seq)
    if reverse:
        scoreList = scoreList[::-1]
    return scoreList

def determine_consensus(name, seq, peaks, qual, median_distance, seed):
    '''
    Aligns and returns the consensus depending on the number of repeats
    If there are multiple peaks, it'll do the normal partial order
    alignment with racon correction
    If there are two repeats, it'll do the special pairwise consensus
    making
    '''
    repeats = ''
    corrected_consensus = ''
    if median_distance > medDistCutoff and len(peaks) >= 1:
        out_F = temp_folder + '/' + name + '_F.fasta'
        out_Fq = temp_folder + '/' + name + '_F.fastq'
        poa_cons = temp_folder + '/' + name + '_consensus.fasta'
        final = temp_folder + '/' + name + '_corrected_consensus.fasta'
        overlap = temp_folder +'/' + name + '_overlaps.sam'
        pairwise = temp_folder + '/' + name + '_prelim_consensus.fasta'
        repeats = split_read(peaks, seq, out_F, qual, out_Fq, name)

        PIR = temp_folder + '/' + name + 'alignment.fasta'
        poa_start = time()
        os.system('%s -read_fasta %s -hb -pir %s \
                  -do_progressive %s &>>poa_messages' \
                  %(poa, out_F, PIR, score_matrix))
        poa_stop = time()
        if timer:
            print('POA took ' + str(poa_stop - poa_start) + ' seconds to run.')
        reads = read_fasta(PIR)

        if repeats == '2':
            Qual_Fasta = open(pairwise, 'w')
            for read in reads:
                if 'CONSENS' not in read:
                    Qual_Fasta.write('>' + read + '\n' + reads[read] + '\n')
            Qual_Fasta.close()
            conspy_start = time()
            os.system('%s %s %s %s >> %s' \
                      %(consensus, pairwise, out_Fq, name, poa_cons))
            conspy_stop = time()
            if timer:
                print('consensus.py took ' + str(conspy_stop - conspy_start) \
                      + ' seconds to run.')

        else:
            for read in reads:
              if 'CONSENS0' in read:
                out_cons_file = open(poa_cons, 'w')
                out_cons_file.write('>' + name + '\n' \
                                    + reads[read].replace('-', '') + '\n')
                out_cons_file.close()

        final = poa_cons
        input_cons = poa_cons
        output_cons = poa_cons.replace('.fasta', '_1.fasta')
        mm_start = time()
        os.system('%s --secondary=no -ax map-ont \
                  %s %s > %s 2> ./minimap2_messages' \
                  %(minimap2, input_cons, out_Fq, overlap))
        mm_stop = time()
        if timer:
            print('minimap2 took ' + str(mm_stop - mm_start) \
                  + ' seconds to run.')
        racon_start = time()
        os.system('%s -q 5 -t 1 \
                  %s %s %s >%s 2>./racon_messages' \
                  %(racon, out_Fq, overlap, input_cons, output_cons))
        racon_stop = time()
        if timer:
            print('Racon took ' + str(racon_stop - racon_start) \
                  + ' seconds to run.')
        final = output_cons
        print(final)
        reads = read_fasta(final)
        for read in reads:
            corrected_consensus = reads[read]
    return corrected_consensus, repeats

def read_fastq_file(seq_file):
    '''
    Takes a FASTQ file and returns a list of tuples
    In each tuple:
        name : str, read ID
        seed : int, first occurrence of the splint
        seq : str, sequence
        qual : str, quality line
        average_quals : float, average quality of that line
        seq_length : int, length of the sequence
    '''
    read_list, lineNum = [], 0
    lastPlus = False
    for line in open(seq_file):
        line = line.rstrip()
        if not line:
            continue
        # make an entry as a list and append the header to that list
        if lineNum % 4 == 0 and line[0] == '@':
            splitLine = line[1:].split('_')
            root, seed = splitLine[0], int(splitLine[1])
            read_list.append([])
            read_list[-1].append(root)
            read_list[-1].append(seed)

        # sequence
        if lineNum % 4 == 1:
            read_list[-1].append(line)

        # quality header
        if lineNum % 4 == 2:
            lastPlus = True

        # quality
        if lineNum % 4 == 3 and lastPlus:
            read_list[-1].append(line)
            avgQ = sum([ord(x)-33 for x in line])/len(line)
            read_list[-1].append(avgQ)
            read_list[-1].append(len(read_list[-1][2]))
            read_list[-1] = tuple(read_list[-1])

        lineNum += 1
    return read_list

def analyze_reads(read_list):
    '''
    Takes reads that are longer than 1000 bases and gives the consensus.
    Writes to R2C2_Consensus.fasta
    '''
    for name, seed, seq, qual, average_quals, seq_length in read_list:
        if seqLenCutoff < seq_length:
            final_consensus = ''
            scoreList = split_SW(name, seed, seq)
            # calculate where peaks are and the median distance between them
            peaks, median_distance = callPeaks(scoreList)

            if figure:
                makeFig(scoreList, peaks, seed, median_distance)

            # make the consensus
            final_consensus, repeats = determine_consensus(name, seq, peaks,
                                                           qual, median_distance,
                                                           seed)
            if final_consensus:
                final_out = open(out_file, 'a')
                final_out.write('>' + name + '_' \
                                + str(round(average_quals, 2)) + '_' \
                                + str(seq_length) + '_' + str(repeats) \
                                + '_' + str(len(final_consensus)))
                final_out.write('\n' + final_consensus + '\n')
                final_out.close()
                os.system('rm ' + temp_folder + '/*')

def main():
    '''Controls the flow of the program'''
    final_out = open(out_file, 'w')
    final_out.close()
    print(input_file)
    read_list = read_fastq_file(input_file)
    analyze_reads(read_list)

if __name__ == '__main__':
    main()
