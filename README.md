# C3POa

[![Github release](https://img.shields.io/github/tag/rvolden/C3POa.svg?label=Version)](https://github.com/rvolden/C3POa/tags)
[![Published in PNAS](https://img.shields.io/badge/Published%20in-PNAS-blue.svg)](https://doi.org/10.1073/pnas.1806447115)
[![GPL license](https://img.shields.io/badge/License-GPL-blue.svg)](http://perso.crans.org/besson/LICENSE.html)

C3POa (**C**oncatemeric **C**onsensus **C**aller with **P**artial **O**rder **a**lignments) is a computational pipeline for calling consensi on R2C2 nanopore data.

**This repository for C3POa is feature frozen. [The most up to date version of C3POa can be found here.](https://github.com/christopher-vollmers/C3POa)**

This version of C3POa changes a lot. The old gonk branch can be found [here](https://github.com/rvolden/C3POa/tree/gonk). The even older version that uses water can be found [here](https://github.com/rvolden/C3POa/tree/water).

## Dependencies

- [Python 3](https://www.python.org/downloads/)
- [NumPy](https://pypi.org/project/numpy/)
- [SciPy](https://pypi.org/project/scipy/)
- [pyabpoa](https://pypi.org/project/pyabpoa/)
- [mappy](https://pypi.org/project/mappy/)
- [tqdm](https://pypi.org/project/tqdm/)
- [Cython](https://pypi.org/project/Cython/)
- [conk](https://github.com/rvolden/conk)
- [racon](https://github.com/isovic/racon)
- [editdistance](https://github.com/roy-ht/editdistance)
- [blat source](https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip) or [blat executable](http://hgdownload.soe.ucsc.edu/admin/exe/)

To fetch and build dependencies, use setup.sh.
setup.sh will download and make the packages that you need to run C3POa (except for blat).
You don't need to have these in your PATH, but if you don't, you'll need to use a [config file](example_config).
The setup script **does not** install programs or add them to your path.
If you use the setup script, you still need to put the paths into a config file (for blat and racon).

```bash
chmod +x setup.sh
./setup.sh
```

Alternatively, you can grab all of the pip installable packages:
```bash
python3 -m pip install --user --upgrade scipy numpy pyabpoa=1.0.5 mappy Cython tqdm setuptools wheel
```
and then build conk and racon manually.

Blat can built from [source](https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip) or you can get an [executable](http://hgdownload.soe.ucsc.edu/admin/exe/).
Please follow the documentation in the blat readme for make instructions.

--------------------------------------------------------------------------------

## Usage

After resolving all of the dependencies, you can run C3POa with python.

## C3POa.py

Preprocessing is now built in.
Preprocessing takes raw 1D nanopore R2C2 reads in fastq (can be zipped) format, removes low quality and short reads and then finds splint sequences in those reads using BLAT.
Preprocessing will also demultiplex reads based on splints that are put into the splint fasta file.
The preprocessor will also look for the alignment psl file in case it was done before.
C3POa won't do the alignment if it finds `output_dir/tmp/splint_to_read_alignments.psl`.
By default, the input file will be chunked into fasta files of `len = group size`.
Use the `-b` option to have chunks be `len = number of reads / number of threads`.

The main algorithmic difference is we only align the splint to the read.
Command line tools have been replaced with their python APIs (except blat and racon).

```bash
python3 C3POa.py -r reads.fastq 
                 -o output/path 
                 -s splint.fasta 
                 -c config 
                 -q 9 
                 -l 1000 
                 -d 500
                 -n 32 
                 -g 1000
```

Arguments:
```
-r  raw reads in fastq format

-o  output path

-s  sequence of DNA splint used in R2C2 protocol in fasta format

-c  config file containing path to BLAT and racon binaries

-q  only reads above this average quality will be retained (9 is recommended)

-l  only reads longer than this number will be retained (1000 recommended)

-d  minimum distance between peaks/minimum insert size. For cDNA, we use 500 (default)

-n  number of threads to use

-g  group size (number of reads given to each thread, default 1000)

-b  split input by number of threads for blat alignment instead of groupSize

-z  use to exclude zero repeat reads

-co compress the output fasta/q files (gzip)

-v  print the C3POa version and exit
```

Example output read (readName_averageQuality_originalReadLength_numberOfRepeats_subreadLength):

```
>efbfbf09-7e2b-48e6-8e57-b3d36886739c_46.53_5798_2_1844
ACAGTCGATCATAGCTTAGCATGCATCGACGATCGATCGATCGA...
```

Example output directory tree:
```
output_dir
├── c3poa.log
├── tmp
│   └── splint_to_read_alignments.psl
├── Splint_1
│   ├── R2C2_Consensus.fasta
│   └── R2C2_Subreads.fastq
└── Splint_2
    ├── R2C2_Consensus.fasta
    └── R2C2_Subreads.fastq
```

--------------------------------------------------------------------------------

## C3POa_postprocessing.py

Trims and reorients consensus sequences generated by C3POa.py to 5'->3' direction.
If given a fasta of oligo dT indexes, it will also demux the reads by index.

```
-i  fasta file containing consensus sequences generated by C3POa.py

-o  directory which output files will be written to

-a  sequence of cDNA adapter sequences in fasta format. Sequence names must be
    3Prime_adapter and 5Prime_adapter

-c  config file containing path to BLAT binary

-x  fasta file of oligo dT indexes

-n  number of threads to use

-u  use to ignore read directionality

-t  use to trim adapters off of the ends of the sequences

-b  use to produce a file with 5' barcode sequences (like 10x reads)

-g  group size (number of reads given to each thread, default 1000)

-bt split input by number of threads instead of groupSize

-co compress the output fasta/q files (gzip)

-v  print the C3POa version and exit
```

```bash
python3 C3POa_postprocessing.py -i /path/to/consensus.fasta -o out_path
                                -c /path/to/config_file -a /path/to/adapter.fasta
```

Example output directory tree with demuxing:
```
output_dir
├── R2C2_oligodT_multiplexing.tsv
├── Index1
│   ├── R2C2_full_length_consensus_reads.fasta
│   ├── R2C2_full_length_consensus_reads_left_splint.fasta
│   └── R2C2_full_length_consensus_reads_right_splint.fasta
├── Index2
│   ├── R2C2_full_length_consensus_reads.fasta
│   ├── R2C2_full_length_consensus_reads_left_splint.fasta
│   └── R2C2_full_length_consensus_reads_right_splint.fasta
└── no_index_found
    ├── R2C2_full_length_consensus_reads.fasta
    ├── R2C2_full_length_consensus_reads_left_splint.fasta
    └── R2C2_full_length_consensus_reads_right_splint.fasta
```
