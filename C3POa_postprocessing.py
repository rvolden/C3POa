#!/usr/bin/env python3
# Roger Volden and Chris Vollmers
# Last updated: 22 May 2018

import sys
import os
import argparse

def argParser():
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

    return vars(parser.parse_args())

args = argParser()
output_path = args['output_path'] + '/'
input_file = args['input_fasta_file']
adapter_file = args['adapter_file']
undirectional = args['undirectional']
trim = args['trim']

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
        # if key not in possible:
        #     raise Exception('Check config file')
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

if args['config'] or args['c']:
    progs = configReader(args['config'])
    blat = progs['blat']
else:
    blat = 'blat'

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

def reverse_complement(sequence):
    '''Returns the reverse complement of a sequence'''
    bases = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', '-':'-'}
    return ''.join([bases[x] for x in list(sequence)])[::-1]

def run_blat(path, infile, adapter_fasta):
    os.system('%s -noHead -stepSize=1 -tileSize=6 -t=DNA q=DNA -minScore=10 \
              -minIdentity=10 -minMatch=1 -oneOff=1 \
              %s %s %s/Adapter_to_consensus_alignment.psl' \
              %(blat,adapter_fasta,infile,path))

def parse_blat(path, reads):
    adapter_dict, iterator = {}, 0

    for name,sequence in reads.items():

        adapter_dict[name] = {}
        adapter_dict[name]['+'] = []
        adapter_dict[name]['-'] = []
        adapter_dict[name]['+'].append(('-', 1, 0))
        adapter_dict[name]['-'].append(('-', 1, len(sequence)))

    for line in open(path + '/Adapter_to_consensus_alignment.psl'):
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

def write_fasta_file(path, adapter_dict, reads):
    out = open(path + 'R2C2_full_length_consensus_reads.fasta', 'w')
    out3 = open(path + 'R2C2_full_length_consensus_reads_left_splint.fasta', 'w')
    out5 = open(path + 'R2C2_full_length_consensus_reads_right_splint.fasta', 'w')
    for name, sequence in reads.items():
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

        if len(plus_list_name) == 1 and len(minus_list_name) == 1:
            if plus_list_position[0] < minus_list_position[0]:
                seq = sequence[plus_list_position[0]:minus_list_position[0]]
                ada = sequence[plus_list_position[0]-40:minus_list_position[0]+40]
                name += '_' + str(len(seq))
                use = False
                if undirectional:
                    direction = '+'
                    use = True
                else:
                    if plus_list_name[0] != minus_list_name[0]:
                        use = True
                        if plus_list_name[0] == '5Prime_adapter':
                            direction = '+'
                        elif plus_list_name[0] == '3Prime_adapter':
                            direction = '-'
                if use:
                    if direction == '+':
                        if trim:
                            out.write('>%s\n%s\n' %(name, seq))
                        else:
                            out.write('>%s\n%s\n' %(name, ada))
                        out5.write('>%s\n%s\n' %(name, reverse_complement(sequence[:plus_list_position[0]])))
                        out3.write('>%s\n%s\n' %(name, sequence[minus_list_position[0]:]))
                    elif direction == '-':
                        if trim:
                            out.write('>%s\n%s\n' %(name, reverse_complement(seq)))
                        else:
                            out.write('>%s\n%s\n' %(name, reverse_complement(ada)))
                        out3.write('>%s\n%s\n' %(name, reverse_complement(sequence[:plus_list_position[0]+40])))
                        out5.write('>%s\n%s\n' %(name, sequence[minus_list_position[0]:]))

def main():
    reads = read_fasta(input_file)

    run_blat(output_path, input_file, adapter_file)
    adapter_dict = parse_blat(output_path, reads)
    write_fasta_file(output_path, adapter_dict, reads)

if __name__ == '__main__':
    main()
