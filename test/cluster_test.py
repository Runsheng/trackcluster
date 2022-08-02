#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/9/2018 10:53 AM
# @Author  : Runsheng     
# @File    : cluster_test.py

from trackcluster.cluster import *
from trackcluster.tracklist import *
import unittest


class ClusterTest(unittest.TestCase):
    def setUp(self):
        # last gene is c24 gene with only single exon
        genes=["unc52", "AT1G06860", "AT2G02100","ATC24-5G67880"]
        gene=genes[3]
        self.gene=gene

        bigg_nano=[]
        with open("./genes/{gene}/{gene}_nano.bed".format(gene=gene)) as f:
            for line_one in f.readlines():
                bigg_one=bigGenePred()
                bigg_one.from_string(line_one)
                bigg_nano.append(bigg_one)

        bigg_gff=[]
        with open("./genes/{gene}/{gene}_gff.bed".format(gene=gene)) as f:
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
        bb=pandas_summary("./genes/unc52/exon_inter.bed")
        print((len(bb)==73093))

    def test_prefilter(self):
        bigg_list_new=prefilter_smallexon(self.bigg_nano, self.bigg_gff, cutoff=50)
        print(( len(self.bigg_nano), len(bigg_list_new) ))

    def test_cal_distance(self):
        D,_=cal_distance(self.bigg, by="ratio")
        print(D)
        D,_=cal_distance(self.bigg, by="ratio_short")
        print(D)
        print(len(D))

    def test_flow(self):
        out="./genes/{gene}/{gene}_simple_coverage.bed".format(gene=self.gene)

        D, bigg_list=flow_cluster(self.bigg_nano, self.bigg_gff, by="ratio_all", intronweight=0.2)
        print(D)
        bigg_nano = add_subread_bigg(bigg_list)

           ### save nessary files
        for bigg in bigg_nano:
            bigg.write_subread()
        bigg_count_write_native(bigg_nano, gff_bed=self.bigg_gff,out=out)


    def test_flow_muti(self):
       D, bigg_list=flow_cluster(self.bigg_nano[1:], self.bigg_gff, by="ratio_all", intronweight=0.2)
       write_D(D, bigg_list, "./genes/unc52/d.csv")

    def tearDown(self):
        self.bigg = None
