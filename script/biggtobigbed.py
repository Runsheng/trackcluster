#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/3/2018 10:55 AM
# @Author  : Runsheng     
# @File    : biggtobigbed.py
"""
todo:Hard-coded path is not used for package!
"""

import os
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-b", "--biggfile",
                    help="the bigg bed file")
parser.add_argument("-o", "--out", default="bigg.bb",
                    help="the output file name")

args = parser.parse_args()

cmd="""
sort -k1,1 -k2,2n {biggs}> {biggs}_s
echo "sort finished"
sed -i '/chrMtDNA/d' {biggs}_s
echo "sed finished"
export PATH="/home/zhaolab1/app/jk/kentUtils/bin/":$PATH
bedToBigBed -as=/home/zhaolab1/myapp/trackcluster/script/bigGenePred.as -type=bed12+8 {biggs}_s /home/zhaolab1/myapp/trackcluster/script/ce10.sizes {out}
""".format(biggs=args.biggfile, out=args.out)

print(cmd)
os.popen(cmd)
