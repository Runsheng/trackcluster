#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/9/2018 10:53 AM
# @Author  : Runsheng     
# @File    : cluster_test.py

from clusterj import *
from tracklist import *
import unittest
from unittest import skip

class ClusterjTest(unittest.TestCase):
    def setUp(self):
        genes=["unc52", "AT1G06860", "AT2G02100", "AT2G43410"]
        gene=genes[3]

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

    def test_get_junction_dic(self):
        site_dic=get_junction_dic(self.bigg)
        print site_dic

    def test_get_start_end_dic(self):
        tss_dic=get_start_end_dic(self.bigg, type="start")
        print tss_dic
        sum_n=0
        for v in tss_dic.values():
            sum_n+=v
        print len(tss_dic), sum_n, len(self.bigg_gff)*5+len(self.bigg_nano)

        tes_dic=get_start_end_dic(self.bigg, type="end")
        print tes_dic

    def test_chain_site(self):
        tes_dic = get_start_end_dic(self.bigg, type="end")
        site_range=chain_site(tes_dic.keys(), offset=10)
        print(site_range, len(tes_dic), len(site_range))

    def test_boundary_correct(self):
        pass

