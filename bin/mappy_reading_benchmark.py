import parse_reads
import mappy as mm

in_file = '../pre/R2C2_Consensus.fasta'
#in_file = '../pre/R2C2_Subreads.fastq.gz'

from time import time
old_start = time()
read_list = parse_reads.read_fasta(in_file)
print('old:', time()-old_start)

mappy_start = time()
read_dict = {}
#new_read_list = []
for read_tuple in mm.fastx_read(in_file, read_comment=False):
    read_dict[read_tuple[0]] = read_tuple[1]
    # read_tuple = list(read_tuple)
    # slen = len(read_tuple[0])
    # read_tuple.append(sum([ord(x)-33 for x in read_tuple[2]])/slen)
    # read_tuple.append(slen)
    # new_read_list.append(tuple(read_tuple))
print('mappy:', time()-mappy_start)
