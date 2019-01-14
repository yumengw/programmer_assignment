# 454 Sequence Process Pipeline

This pipeline is designed to remove primer and adaptor sequences based on blast. The example input fasta and quality files are from a dataset generated using 454 sequencing platform. The script will blast the dataset against the given primer and adaptor sequences and generate output in m8 format.

## Prerequsite

*python2.7 or python3.6 (or above)

*(Biopython)[https://biopython.org]
To install Biopython through pip:
```
pip install biopython
```
*(blast command line)[https://blast.ncbi.nlm.nih.gov/Blast.cgi]
To install blast, follow the instruction (here)[https://www.ncbi.nlm.nih.gov/books/NBK279690/]

## Required input
Fasta and quality files containing the reads and corresponding quality scores.

## Commands for running
```
python SeqProcessor.py -r Files_for_test/test.fna -q Files_for_test/test.qual
```

## Help
use *python SeqProcessor.py* to look for all options

## Output
Adaptor and primer sequences are:
 
Primer Sequence:CGCCGTTTCCCAGTAGGTCTC
Adaptor Sequence:ACTGAGTGGGAGGCAAGGCACACAGGGGATAGG
 
The program should generate the following output:

1) Total number of reads in the dataset.
2) Total number of reads greater than 100 bp.
3) Total number of reads with average quality scores greater than 20.
4) Total number of reads with primer sequences.
5) Total number of reads with adaptor sequences.
6) Total number of reads with both primer and adaptor sequences.

In addition, your program needs to generate the following files:

1) Blast output file in m8 format.
2) Fasta file containing reads greater than 100bp, average read quality scores greater than 20, primers and adaptors trimmed.
3) Tab de-limited text file containing the read identifiers along with the starting and end positions of the primer or adaptor sequences.

