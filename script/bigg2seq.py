#!/usr/bin/env python
#-*- coding: utf-8 -*-
# @Time    : 19/2/2024 2:14pm
# @Author  : Runsheng
# @File    : bigg2seq.py

"""
This script is used to write the exon information from the bigg file
can be used to output the
transcript sequence, CDS sequence and ensemble like EXON(uppercase)+intron(lowercase or N) sequence.

"""
import argparse
import os,sys,inspect
from trackcluster.utils import fasta2dic

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(os.path.dirname(currentdir))
sys.path.insert(0,parentdir)

from trackcluster.tracklist import read_bigg, write_bigg

parser=argparse.ArgumentParser()
parser.add_argument("-b", "--biggfile",
                    help="the bigg bed file")
parser.add_argument("-r", "--reference",
                    help="the genome reference")
parser.add_argument("-o", "--out", default="bigg.fasta",
                    help="the output file name, default is bigg.fasta")
parser.add_argument("-m", "--mode", default="exon",
                    help="the format of the output, can be 'exon', 'cds', "
                         "or 'ensembl': ensemble EXON(uppercase)+intron(lowercase) sequence"
                         "default is 'exon' mode")

args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

# make a file using the functions
outfile=args.out
refdic=fasta2dic(args.reference)

bigg_l=read_bigg(args.biggfile)
with open (outfile, "w") as fw:
    if args.mode=="exon":
        for bigg in bigg_l:
            bigg.get_exon()
            bigg.bind_chroseq(refdic, gap=0, intron=False)
            name=bigg.name
            seq=bigg.seq_chro
    elif args.mode=="cds":
        for bigg in bigg_l:
            bigg.get_exon()
            bigg.get_cds()
            name=bigg.name
            seq=bigg.seq_cds

    for bigg in bigg_l:
        bigg.get_exon()
        fw.write(bigg.name+"\t"+str(bigg.exonlen)+"\t"+bigg.geneName+"\t"+bigg.ttype+"\n")