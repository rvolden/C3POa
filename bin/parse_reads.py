#!/usr/bin/env python3
# Roger Volden

def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    if inFile.endswith('.gz'):
        import gzip
        fh = gzip.open(inFile, 'rt')
    else:
        fh = open(inFile)
    readDict = {}
    for line in fh:
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            if readDict:
                readDict[lastHead] = ''.join(readDict[lastHead])
            readDict[line[1:]] = []
            lastHead = line[1:]
        else:
            readDict[lastHead].append(line)
    if readDict:
        readDict[lastHead] = ''.join(readDict[lastHead])
    return readDict

def read_fastq(inFile):
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
    if inFile.endswith('.gz'):
        import gzip
        fh = gzip.open(inFile, 'rt')
    else:
        fh = open(inFile)

    read_list, lineNum = [], 0
    for line in fh:
        line = line.rstrip()
        if not line:
            continue
        # make an entry as a list and append the header to that list
        if lineNum % 4 == 0 and line[0] == '@':
            splitLine = line[1:].split('_')
            # Kayla: edited to handle hairpin split reads with pre and post
            if len(splitLine) == 2:
                root, seed = splitLine[0], int(splitLine[1])
            else:
                root, seed = splitLine[0] + '_' + splitLine[1], int(splitLine[2])
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
            lastPlus = False

        lineNum += 1

    return read_list

