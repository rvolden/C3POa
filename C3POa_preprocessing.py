import os
import numpy as np
import argparse


parser=argparse.ArgumentParser()
parser.add_argument('-i','--input_fastq_file',type=str)
parser.add_argument('-o','--output_path',type=str)
parser.add_argument('-q','--quality_cutoff',type=str)
parser.add_argument('-l','--read_length_cutoff',type=str)
parser.add_argument('-s','--splint_file',type=str)


args=parser.parse_args()
output_path=args.output_path+'/'	         
input_file = args.input_fastq_file
quality_cutoff = float(args.quality_cutoff)
read_length_cutoff = float(args.read_length_cutoff)
splint_file = args.splint_file

blat='blat'

def read_and_filter_fastq(input_file):
    length=0
    for line in open(input_file):
        length+=1
    iterator=0
    read_dict={}
    input1=open(input_file,'r')
    while iterator<length:
        a=input1.readline().strip()
        b=input1.readline().strip()
        c=input1.readline().strip()
        d=input1.readline().strip()
        quals=[]
        for character in d:
            number = ord(character)-33
            quals.append(number) 
        average=np.average(quals)


        name=a[1:].split()[0]
        sequence=b
        qual=d

        if len(sequence)>=read_length_cutoff and average>=quality_cutoff: 
            read_dict[name]=(sequence,qual)
        iterator+=4

    return(read_dict)

def run_blat(path,reads):

    fasta_file=path+'/R2C2_temp_for_BLAT.fasta'
    input_file_fasta=open(fasta_file,'w')

    for read in reads:
        input_file_fasta.write('>'+read+'\n'+reads[read][0]+'\n')
    input_file_fasta.close()

    os.system('%s -noHead -stepSize=1 -t=DNA q=DNA -minScore=15 -minIdentity=10 %s %s %s/Splint_to_read_alignments.psl' %(blat,splint_file,fasta_file,path))

def parse_blat(path):
    alignment_file=path+'/Splint_to_read_alignments.psl'
    fasta_file=path+'/R2C2_temp_for_BLAT.fasta'
    adapter_dict={}

    length=0
    for line in open(fasta_file):
        length+=1

    iterator=0
    infile=open(fasta_file,'r')
    while iterator<length:
        line=infile.readline()
        sequence=infile.readline()
        name=line[1:].strip()
        adapter_dict[name]={}

        adapter_dict[name]['+']=[]
        adapter_dict[name]['-']=[]

        adapter_dict[name]['+'].append(('-',1,0))
        adapter_dict[name]['-'].append(('-',1,len(sequence)))

        iterator+=2   


    burn_dict={}
    for line in open(alignment_file):
        a=line.strip().split('\t')
        read_name=a[9]
        adapter=a[13]
        strand=a[8]
        gaps=float(a[5])
        score=float(a[0])
        sequence_length=float(a[10])  
        if gaps<50 and score>50:   # Looks for unspliced quality alignment
          if strand=='+':
                start=int(a[11])-int(a[15])
                end=int(a[12])+(int(a[14])-int(a[16]))
          if strand=='-':
                start=int(a[11])-(int(a[14])-int(a[16]))
                end=int(a[12])+int(a[15])

          position=min(max(0,int(start+((end-start)/2))),sequence_length-1)

          adapter_dict[read_name][strand].append((adapter,float(a[0]),position))


    return adapter_dict


def write_fastq_files(path,adapter_dict,reads):

    success=0
    os.system('mkdir '+path+'/splint_reads')
    for read in reads:
        name=read
        sequence=reads[read][0]
        quality=reads[read][1]

        adapter_plus=sorted(adapter_dict[name]['+'],key=lambda x: x[2],reverse=False)
        adapter_minus=sorted(adapter_dict[name]['-'],key=lambda x: x[2],reverse=False)

        plus_list_name=[]
        plus_list_position=[]
        minus_list_name=[]
        minus_list_position=[]


        for adapter in adapter_plus:
            if adapter[0]!='-':
                plus_list_name.append(adapter[0])
                plus_list_position.append(adapter[2])
        for adapter in adapter_minus:
            if adapter[0]!='-':
                minus_list_name.append(adapter[0])
                minus_list_position.append(adapter[2])

        if 'Splint' in plus_list_name or 'Splint' in minus_list_name:
            success+=1
            try:
                out_fastq=open(path+'/splint_reads/'+str(int(success/4000))+'/R2C2_raw_reads'+'.fastq','a')
            except:
                os.system('mkdir '+path+'/splint_reads/'+str(int(success/4000)))
                out_fastq=open(path+'/splint_reads/'+str(int(success/4000))+'/R2C2_raw_reads'+'.fastq','w')

        
            if len(plus_list_name)>0:
                out_fastq.write('>'+name+'_'+str(plus_list_position[0])+'\n'+sequence+'\n+\n'+quality+'\n')
            else:
                out_fastq.write('>'+name+'_'+str(minus_list_position[0])+'\n'+sequence+'\n+\n'+quality+'\n')



def main():
    print('Reading and filtering fastq file')
    reads=read_and_filter_fastq(input_file)
    print('Running BLAT to find splint locations (This can take hours)')
    run_blat(output_path,reads)
    print('Parsing BLAT output')
    adapter_dict=parse_blat(output_path)
    print('Writing fastq output files in bins of 4000 into separate folders')
    write_fastq_files(output_path,adapter_dict,reads)

main()
