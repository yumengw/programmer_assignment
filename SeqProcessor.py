#!/usr/bin/env python

import time
import sys
import os
import numpy as np
from optparse import OptionParser
from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqRecord import SeqRecord

parser = OptionParser()
parser.add_option("-r", "--read", dest="read_file",
                  help="raw reads fasta file", metavar="FILE")
parser.add_option("-q", "--qual", dest="qual_file",
                  help="raw reads quality file", metavar="FILE")
parser.add_option("-l", "--len", dest="length_cutoff",
                  help="length cutoff", type="int", default=100)
parser.add_option("-s", "--score", dest="qual_cutoff",
                  help="length cutoff", type="int", default=20)
parser.add_option("-p", "--primer", dest="primer_seq",
                  help="primer sequence", type="string", default="CGCCGTTTCCCAGTAGGTCTC")
parser.add_option("-a", "--adaptor", dest="adaptor_seq",
                  help="adaptor sequence", type="string", default="ACTGAGTGGGAGGCAAGGCACACAGGGGATAGG")
parser.add_option("-w", "--word_size", dest="word_size",
                  help="blast word size", type="int", default=7)
parser.add_option("-e", "--evalue", dest="evalue",
                  help="blast evalue", type="float", default=0.001)
parser.add_option("-t", "--tolerate", dest="tolerate",
                  help="maximum tolerate error rate when triming primer and adaptor", type="float", default=0.1)
parser.add_option("-o", "--outdir", dest="outdir", metavar="STRING",
                  help="output dir", type="string", default="Results")

class SeqProcessor():
    def __init__(self):
        '''
        Initialize from input options, and initialize count for qualified reads
        '''
        print('Initilizing input options...')
        self.seq_file = options.read_file
        self.qual_file = options.qual_file
        self.primer_seq = options.primer_seq
        self.adaptor_seq = options.adaptor_seq
        self.blast_out = None

        self.n_reads_total = 0
        self.n_reads_gt_100 = 0
        self.n_reads_avg_qual_gt_20 = 0
        self.n_reads_with_primer = 0
        self.n_reads_with_adaptor = 0
        self.n_reads_with_primer_adaptor = 0
        
    def blast_primer_adaptor(self):
        '''
        Run blastn using given parameters, and update self.blast_out
        Usage: SeqProcessor.blast_primer_adaptor()
        '''
        query_string  = '>primer' + '\n' + str(self.primer_seq) + '\n' + \
                        '>adaptor' + '\n' + str(self.adaptor_seq)
        if not os.path.exists(options.outdir): # generate outdir if not exist
            os.makedirs(options.outdir)
        num_seq = len([1 for line in open(self.seq_file) if line.startswith(">")])
        blast_cmd = NcbiblastnCommandline(subject = self.seq_file, word_size = options.word_size, max_target_seqs = num_seq*2, \
                                          evalue = options.evalue, outfmt=6, out=options.outdir+os.sep+'blast_out.m8')
        print('Runing blast against primer and adaptor...')
        blast_cmd(stdin=query_string)
        # update self.blast_out as a dictionary {seq_id: alignment}
        self.blast_out = {'primer':{}, 'adaptor':{}}
        try: blast_results = open(options.outdir+os.sep+'blast_out.m8', 'r')
        except:
            print('Error runing blast, please check your input and make sure you installed Blast locally.')
            sys.exit(1)
        for line in blast_results:
            item = line.rstrip().split('\t')
            seq_id = item[1]
            query_id = item[0]
            if query_id=='primer':
                try: self.blast_out['primer'][seq_id].append(item)
                except: self.blast_out['primer'][seq_id] = [item]
            elif query_id=='adaptor':
                try: self.blast_out['adaptor'][seq_id].append(item)
                except: self.blast_out['adaptor'][seq_id] = [item]
                
    def count_qualified_seq(self):
        '''
        Update count information based on quality file and blast result
        Usage: SeqProcessor.count_qualified_seq()
        '''
        if self.blast_out == None:
            print('Please run blast_primer_adaptor() before filtering sequences')
            sys.exit(1)
        try: PairedFastaQualIterator(open(self.seq_file), open(self.qual_file))
        except:
            print('Error reading sequence file and matched quality file, please double check')
            sys.exit(1)
            
        outfile_filtered_seq = open(options.outdir+os.sep+'filtered_seq.fna', 'w')
        out_align = open(options.outdir+os.sep+'alignment.tsv', 'w')
        out_align.write('seq_id\tstart\tend\talign_to\n')
        for record in PairedFastaQualIterator(open(self.seq_file), open(self.qual_file)):
            self.n_reads_total += 1
            if self.n_reads_total % 1000 == 0: print("processing read %d ..." %(self.n_reads_total))
            ## search and trim primer and adptor
            has_primer, primer_start, primer_end = self._trim_primer(record)
            has_adaptor, adaptor_start, adaptor_end = self._trim_adaptor(record)
            trimed_start = 0
            trimed_end = len(record.seq)
            if has_primer:
                self.n_reads_with_primer += 1
                trimed_start = primer_end
                out_align.write(record.id +'\t'+ str(primer_start) +'\t'+ str(primer_end) +'\t'+ 'primer\n')
                
            if has_adaptor:
                self.n_reads_with_adaptor += 1
                trimed_end = adaptor_start - 1
                out_align.write(record.id +'\t'+ str(adaptor_start) +'\t'+ str(adaptor_end) +'\t'+ 'adaptor\n')
                
            if has_primer and has_adaptor:
                self.n_reads_with_primer_adaptor += 1
                
            trimed_seq = record[trimed_start:trimed_end]
            if len(trimed_seq.seq) > options.length_cutoff: self.n_reads_gt_100 += 1
            if np.mean(trimed_seq.letter_annotations["phred_quality"]) > options.qual_cutoff:
                self.n_reads_avg_qual_gt_20 += 1
            if len(trimed_seq.seq) > options.length_cutoff and \
                np.mean(trimed_seq.letter_annotations["phred_quality"]) > options.qual_cutoff:
                SeqIO.write(trimed_seq, outfile_filtered_seq, "fasta")
    
    def _trim_primer(self, record):
        if record.id in self.blast_out['primer']:
            for each_align in self.blast_out['primer'][record.id]:
                max_error = int(round(int(each_align[3])*options.tolerate))
                if int(each_align[8]) > max_error: continue # not 5' anchored
                if int(each_align[4]) + int(each_align[5]) > max_error: continue
                return 1, int(each_align[8]), int(each_align[9])
        return 0, -1, -1
    
    def _trim_adaptor(self, record):
        if record.id in self.blast_out['adaptor']:
            for each_align in self.blast_out['adaptor'][record.id]:
                max_error = int(round(int(each_align[3])*options.tolerate))
                if int(each_align[4]) + int(each_align[5]) > max_error: continue
                return 1, int(each_align[8]), int(each_align[9])
        return 0, -1, -1
    
    def summary(self):
        print("------------ Summary -------------")
        print("Total number of reads in the dataset: %d" %self.n_reads_total)
        print("Total number of reads greater than %d bp: %d" %(options.length_cutoff, self.n_reads_gt_100))
        print("Total number of reads with average quality scores greater than %d: %d" \
                %(options.qual_cutoff, self.n_reads_avg_qual_gt_20))
        print("Total number of reads with primer sequences: %d" %self.n_reads_with_primer)
        print("Total number of reads with adaptor sequences: %d" %self.n_reads_with_adaptor)
        print("Total number of reads with both primer and adaptor sequences: %d" %self.n_reads_with_primer_adaptor)
    
if __name__ == "__main__":

    ## check python version
    if sys.version < "2.7":
        raise Exception("Requires Python 2.7 or later.")

    ## load options
    (options, args) = parser.parse_args()
    if len(sys.argv) < 2:
        parser.print_help()
        exit(0)

    ## start calculation
    start_time = time.time()
    ## init processor based on input options
    seq_processor = SeqProcessor()
    ## run blast
    seq_processor.blast_primer_adaptor()
    ## based on blast, trim adaptor and primer
    seq_processor.count_qualified_seq()
    print("\n--- Job finished in %s seconds ---\n" % (time.time() - start_time))
    ## output summary
    seq_processor.summary()
