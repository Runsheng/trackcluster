#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os,sys,inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(os.path.dirname(currentdir))
sys.path.insert(0,parentdir)

from trackcluster.convert import gff_to_bigGenePred
from trackcluster.gff import GFF
from trackcluster.tracklist import write_bigg

parser=argparse.ArgumentParser()
parser.add_argument("-i", "--gff",
                    help="the sorted gff file")
parser.add_argument("-o", "--out", default="bigg.bed",
                    help="the output bigGenePred file name")

parser.add_argument("-k", "--key", default="ID",
                    help="The key used as gene name in gff line attr column, default is ID, which is used in most Ensembl gff."
                         "like ID=gene:ENSMUSG00000064842;Name=let-1")

args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

# make a file using the functions

gff = GFF(args.gff)
bigg_list = gff_to_bigGenePred(gff, indicator=args.key)

write_bigg(bigg_list, args.out)

