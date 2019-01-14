# 454 Sequence Process Pipeline

This pipeline is designed to remove primer and adaptor sequences based on blast. The example input fasta and quality files are from a dataset generated using 454 sequencing platform. The script will blast the dataset against the given primer and adaptor sequences and generate output in m8 format.

## Prerequsite

1. python2.7 or python3.6 (or above)

2. [Biopython](https://biopython.org)
To install Biopython through pip:
```
pip install biopython
```

3. [blast command line](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
To install blast, follow the instruction [here](https://www.ncbi.nlm.nih.gov/books/NBK279690/)

## Required input
Fasta and quality files containing the reads and corresponding quality scores.
Default Primer Sequence:CGCCGTTTCCCAGTAGGTCTC
Default Adaptor Sequence:ACTGAGTGGGAGGCAAGGCACACAGGGGATAGG

## Commands for running
```
python SeqProcessor.py -r Files_for_test/test.fna -q Files_for_test/test.qual
```

## Help
```
python SeqProcessor.py

Usage: SeqProcessor.py [options]

Options:
  -h, --help            show this help message and exit
  -r FILE, --read=FILE  raw reads fasta file
  -q FILE, --qual=FILE  raw reads quality file
  -l LENGTH_CUTOFF, --len=LENGTH_CUTOFF
                        length cutoff
  -s QUAL_CUTOFF, --score=QUAL_CUTOFF
                        length cutoff
  -p PRIMER_SEQ, --primer=PRIMER_SEQ
                        primer sequence
  -a ADAPTOR_SEQ, --adaptor=ADAPTOR_SEQ
                        adaptor sequence
  -w WORD_SIZE, --word_size=WORD_SIZE
                        blast word size
  -e EVALUE, --evalue=EVALUE
                        blast evalue
  -t TOLERATE, --tolerate=TOLERATE
                        maximum tolerate error rate when triming primer and
                        adaptor
  -o STRING, --outdir=STRING
                        output dir
```

## Output

The program should generate the following output:

1) Total number of reads in the dataset.
2) Total number of reads greater than 100 bp (cutoff can be changed through -l).
3) Total number of reads with average quality scores greater than 20 (cutoff can be changed through -s).
4) Total number of reads with primer sequences.
5) Total number of reads with adaptor sequences.
6) Total number of reads with both primer and adaptor sequences.

In addition, the script will generate a "Results" folder (folder name can be change through -o) containing the following files:

1) blast_out.m8: Blast output file in m8 format.
2) filtered_seq.fna: Fasta file containing reads greater than 100bp, average read quality scores greater than 20, primers and adaptors trimmed.
3) alignment.tsv: Tab de-limited text file containing the read identifiers along with the starting and end positions of the primer or adaptor sequences.

