#!/usr/bin/env python
#-*- coding: utf-8 -*-

from pysam import AlignmentFile
import argparse
import os,sys,inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(os.path.dirname(currentdir))
sys.path.insert(0,parentdir)

# self import
from trackcluster.flow import flow_bamconvert

parser=argparse.ArgumentParser()
parser.add_argument("-b", "--bamfile",
                    help="the sorted and indexed bam file")
parser.add_argument("-o", "--out", default="bigg.bed",
                    help="the output file name")
parser.add_argument("-s", "--score",  type=int, default=30,
                    help="The min mapq score used to keep a read")

parser.add_argument("-g", "--group", default=None,
                    help="The name used for the group of this sample, if not given, use the prefix of the bamfile")
parser.add_argument("-d", "--wkdir", default=None,
                    help="The dir path contain the file, if not given, use the current dir")

args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

# default handler
wkdir=os.getcwd() if args.wkdir is None else args.wkdir
prefix=args.bamfile.split(".")[0] if args.group is None else args.group

flow_bamconvert(wkdir=wkdir,bamfile=args.bamfile,out=args.out, prefix=prefix, score=int(args.score))
