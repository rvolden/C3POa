import sys
import os
import numpy as np
import argparse
import editdistance as ld

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_fasta_file', type=str)
parser.add_argument('-o', '--output_path', type=str)
parser.add_argument('-n', '--nextera_index_file', type=str)
parser.add_argument('-t', '--tso_index_file', type=str)

args = parser.parse_args()
output_path = args.output_path + '/'
input_file = args.input_fasta_file
nextera_file = args.nextera_index_file
tso_file = args.tso_index_file

def read_fasta(inFile):
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

def reverse_complement(sequence):
    bases = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', '-':'-'}
    return ''.join([bases[x] for x in list(sequence)])[::-1]

def demultiplex(reads, Nextera_Indexes, TSO_Indexes):
    Nexts = read_fasta(Nextera_Indexes)
    TSOs = read_fasta(TSO_Indexes)
    indexed_reads = {}
    counter = 0
    for read, complete_sequence in reads.items():
        if counter % 10000 == 0:
            print(str(counter) + ' of ' + str(len(reads)))
        if len(complete_sequence) > 300:
                sequence = complete_sequence[:300]

                TSO_list, Nextera_list = [], []
                max_dist1, Next_match = 4, ''
                for Next in Nexts:
                    hseq = Nexts[Next]
                    length = len(hseq)
                    temp_dists = []
                    for i in range(len(sequence) - length):
                        dist = ld.eval(hseq, sequence[i:i+length])
                        temp_dists.append(dist)
                    smallest_dist = sorted(temp_dists)[0]
                    Nextera_list.append((Next, smallest_dist))

                max_dist2, TSO_match = 4, ''
                for TSO in TSOs:
                    hseq = TSOs[TSO]
                    length = len(hseq)
                    temp_dists = []
                    for i in range(len(sequence) - length):
                        dist = ld.eval(hseq, sequence[i:i+length])
                        temp_dists.append(dist)
                    smallest_dist = sorted(temp_dists)[0]
                    TSO_list.append((TSO, smallest_dist))

                sorted_Nextera_list = sorted(Nextera_list, key=lambda x: x[1])
                sorted_TSO_list = sorted(TSO_list, key=lambda x: x[1])
                if sorted_Nextera_list[0][1] < max_dist1 and \
                   sorted_Nextera_list[0][1] < sorted_Nextera_list[1][1]-1:
                    Next_match = sorted_Nextera_list[0][0]
                if sorted_TSO_list[0][1] < max_dist2 and \
                   sorted_TSO_list[0][1] < sorted_TSO_list[1][1]-1:
                    TSO_match = sorted_TSO_list[0][0]

                read += '|' + Next_match + '_' + TSO_match
                indexed_reads[read] = complete_sequence
                counter += 1
    return indexed_reads

def write_fasta_file(path, reads):
    out = open(path + '/Indexed_reads.fasta', 'w')
    for name, sequence in reads.items():
        out.write('>%s\n%s\n' %(name, sequence))

def main():
    reads = read_fasta(input_file)
    indexed_reads = demultiplex(reads, nextera_file, tso_file)
    write_fasta_file(output_path, indexed_reads)

main()
