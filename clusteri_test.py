#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/9/2018 10:53 AM
# @Author  : Runsheng     
# @File    : cluster_test.py

from clusteri import *
import unittest
from unittest import skip

class ClusterTest(unittest.TestCase):
    def setUp(self):
        bigg_nano=[]
        with open("./test/unc52_sw.bed") as f:
            for line_one in f.readlines():
                bigg_one=bigGenePred()
                bigg_one.from_string(line_one)
                bigg_nano.append(bigg_one)

        bigg_gff=[]
        with open("./test/unc52_gff.bed") as f:
            for line_one in f.readlines():
                bigg_one=bigGenePred()
                bigg_one.from_string(line_one)
                bigg_gff.append(bigg_one)
        self.bigg_nano=bigg_nano
        self.bigg_gff=bigg_gff
        self.bigg=bigg_nano+bigg_gff

        for i in self.bigg:
            i.to_bedfile(rewrite=False)


    def test_IO(self):
        pass
        #self.assertEquals(len(self.bigg), 323)

    def test_prefilter(self):
        bigg_list_new=prefilter_smallexon(self.bigg_nano, self.bigg_gff, cutoff=50, core=40)
        print len(self.bigg_nano), len(bigg_list_new)


    def test_muti_wrapper(self):
        pair_list=[]

        for i in self.bigg_nano:
            pair_list.append((i, self.bigg_gff))

        #print wrapper_bedtools_intersection_muti(pair_list, use="exon", core=40)

    def test_cluster(self):
        sample = self.bigg[0:]
        D,bigg_list=cal_distance(sample,intronweight=0.5,by="ratio", core=40)
        print(D)
        #keep=filter_D(D, sample)
        #print len(bigg_list)


    def tearDown(self):
        self.bigg = None
