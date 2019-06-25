#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/9/2018 10:53 AM
# @Author  : Runsheng     
# @File    : cluster_test.py

from cluster import *
from tracklist import *
import unittest
from unittest import skip

class ClusterTest(unittest.TestCase):
    def setUp(self):
        genes=["unc52", "AT1G06860", "AT2G02100"]
        gene=genes[1]

        bigg_nano=[]
        with open("./test/genes/{gene}/{gene}_nano.bed".format(gene=gene)) as f:
            for line_one in f.readlines():
                bigg_one=bigGenePred()
                bigg_one.from_string(line_one)
                bigg_nano.append(bigg_one)

        bigg_gff=[]
        with open("./test/genes/{gene}/{gene}_gff.bed".format(gene=gene)) as f:
            for line_one in f.readlines():
                bigg_one=bigGenePred()
                bigg_one.from_string(line_one)
                bigg_gff.append(bigg_one)

        self.bigg_nano=bigg_nano
        self.bigg_gff=bigg_gff
        self.bigg=bigg_nano+bigg_gff

    def test_IO(self):
        pass
        #self.assertEquals(len(self.bigg), 323)

    def test_pandas_summary(self):
        bb=pandas_summary("./test/genes/unc52/exon_inter.bed")
        print(len(bb)==73093)

    def test_prefilter(self):
        bigg_list_new=prefilter_smallexon(self.bigg_nano, self.bigg_gff, cutoff=50)
        print len(self.bigg_nano), len(bigg_list_new)

    def test_cal_distance(self):
        D,_=cal_distance(self.bigg)
        print D

    def test_flow(self):
       D, bigg_list=flow_cluster(self.bigg_nano, self.bigg_gff, by="ratio_all", intronweight=0.5)
       bigg_nano = add_subread_bigg(bigg_list)

       ### save nessary files
       for bigg in bigg_nano:
           bigg.write_subread()

       bigg_count_write(bigg_nano, out="./test/unc_52_simple_coverage.bed")


    def test_flow_muti(self):
       D, bigg_list=flow_cluster(self.bigg_nano[1:100], self.bigg_gff, by="ratio_all", intronweight=0.2)
       write_D(D, bigg_list, "./test/d.csv")

    def tearDown(self):
        self.bigg = None
