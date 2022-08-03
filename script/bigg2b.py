#!/usr/bin/env python

#"""
#"""

import os
import argparse
import sys,inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(os.path.dirname(currentdir))
sys.path.insert(0,parentdir)


parser=argparse.ArgumentParser()
parser.add_argument("-b", "--biggfile",
                    help="the bigg bed file")
parser.add_argument("-o", "--out", default="bigg.bb",
                    help="the output file name")
parser.add_argument("-k", "--kent", default="/home/bin",
                    help="the kent bin location")
parser.add_argument("-s", "--sizefile", default=currentdir+"/ce10.sizes",
                    help="the genome size file")

args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

cmd="""
LC_COLLATE=C sort -k1,1 -k2,2n {biggs}> {biggs}_s
echo "sort finished"
export PATH={kent}:$PATH
bedToBigBed -as={currentdir}/bigGenePred.as -type=bed12+8 {biggs}_s {sizefile} {out}
""".format(biggs=args.biggfile, out=args.out, currentdir=currentdir, sizefile=args.sizefile,
           kent=args.kent)

print(cmd)
os.popen(cmd)
