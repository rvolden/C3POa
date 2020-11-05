#!/usr/bin/env python3
# Roger Volden and Chris Vollmers

import sys
import os
import argparse
import mappy as mm

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

def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readDict = {}
    for read in mm.fastx_read(inFile, read_comment=False):
        readDict[read[0]] = read[1]
    return readDict

def run_blat(path, infile, adapter_fasta, blat):
    align_psl = path + 'adapter_to_consensus_alignment.psl'
    os.system('{blat} -noHead -stepSize=1 -tileSize=6 -t=DNA -q=DNA -minScore=10 \
              -minIdentity=10 -minMatch=1 -oneOff=1 {adapters} {reads} {psl}'\
              .format(blat=blat, adapters=adapter_fasta, reads=infile, psl=align_psl))

def parse_blat(path, reads):
    adapter_dict, iterator = {}, 0

    for name, sequence in reads.items():
        adapter_dict[name] = {}
        adapter_dict[name]['+'] = []
        adapter_dict[name]['-'] = []
        adapter_dict[name]['+'].append(('-', 1, 0))
        adapter_dict[name]['-'].append(('-', 1, len(sequence)))

    for line in open(path + 'adapter_to_consensus_alignment.psl'):
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

def write_fasta_file(args, adapter_dict, reads):
    path = args.output_path
    undirectional = args.undirectional
    barcoded = args.barcoded

    out = open(path + 'R2C2_full_length_consensus_reads.fasta', 'w')
    out3 = open(path + 'R2C2_full_length_consensus_reads_left_splint.fasta', 'w')
    out5 = open(path + 'R2C2_full_length_consensus_reads_right_splint.fasta', 'w')
    if barcoded:
        out10X = open(path + 'R2C2_full_length_consensus_reads_10X_sequences.fasta', 'w')

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

    reads = read_fasta(args.input_fasta_file)

    run_blat(args.output_path, args.input_fasta_file, args.adapter_file, blat)
    adapter_dict = parse_blat(args.output_path, reads)
    write_fasta_file(args, adapter_dict, reads)

if __name__ == '__main__':
    args = parse_args()
    main(args)
