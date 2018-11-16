#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/1/2018 12:03 AM
# @Author  : Runsheng     
# @File    : track_test.py


import unittest

from track import *
from utils import fasta2dic

class PlotTest(unittest.TestCase):
    def setUp(self):
        bigg=[]
        with open("./test/unc52_sw.bed") as f:
            for line_one in f.readlines():
                bigg_one=bigGenePred()
                bigg_one.from_string(line_one)
                bigg.append(bigg_one)
        with open("./test/unc52_gff.bed") as f:
            for line_one in f.readlines():
                bigg_one=bigGenePred()
                bigg_one.from_string(line_one)
                bigg.append(bigg_one)
        self.bigg=bigg

    def test_IO(self):
        sample=self.bigg[0]
        sample.to_bedstr()
        print sample.exon_str

    def test_cal_distance(self):

        bed1, bed2=self.bigg[0:2]
        #print(bed1.bedfile_cal_distance_exon(bed2))
        #print(bed1.bedfile_cal_distance_intron(bed2))
        pass


    def test_bindseq(self):
        ref_dict=fasta2dic("/home/zhaolab1/reference/ce10.fa")
        bigg_one=self.bigg[-1]

        bigg_one.bind_chroseq(ref_dict, gap=100, intron=True)
        print bigg_one.seq_chro
        print bigg_one


    def __test_orfs(self):
        """
        can only run with ce10 ref
        """
        ref_dict=fasta2dic("/home/zhaolab1/reference/ce10.fa")
        bigg_one=self.bigg[30]
        bigg_one.bind_chroseq(ref_dict, gap=0, intron=False)
        print bigg_one.seq_chro
        ans=bigg_one.find_orfs_with_trans()

        print ans
        print bigg_one

    def test_orf(self):
        pass





