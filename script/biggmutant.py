#!/usr/bin/env python
#-*- coding: utf-8 -*-

import argparse
import os,sys,inspect

from trackcluster.tracklist import bigg_addvalue

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(os.path.dirname(currentdir))
sys.path.insert(0,parentdir)


"""
the script used to add a simple value to a bigglist,
like score, GeneName, GeneName2(group), ttype(isoform_anno, read_nano, region_mark)
pos is designed to be score, GeneName, GeneName2(group), ttype(isoform_anno, read_nano, region_mark), name2 (subreads)
"""

parser=argparse.ArgumentParser()
parser.add_argument("-i", "--input",
                    help="the bigg input file")
parser.add_argument("-o", "--out", default="biggnew.bed",
                    help="the output file name")
parser.add_argument("-k", "--key",  type=int, default=30,
                    help="The position you want to change,  can only be one of score, genename, group, type, subread")
parser.add_argument("-v", "--value", default=None,
                    help="The new value you want to give to the key position")

args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

name_dic={"score":"score",
          "genename": "GeneName",
          "group":"GeneName2",
          "type":"ttype",
          "subread":"name2"}

try:
    key=name_dic[args.key]
except KeyError:
    print("key must be in one of score, genename, group, type, subread")

bigg_addvalue(bigg_file=args.input, value_str=args.value, poskey=args.key, out=args.out)