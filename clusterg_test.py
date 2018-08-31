#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/9/2018 10:53 AM
# @Author  : Runsheng     
# @File    : cluster_test.py

from clusterg import *
import unittest
from unittest import skip

class ClusterTest(unittest.TestCase):
    def setUp(self):
        bigg = []
        with open("./test/unc52_sw.bed") as f:
            for line_one in f.readlines():
                bigg_one = bigGenePred()
                bigg_one.from_string(line_one)
                bigg.append(bigg_one)
        self.bigg = bigg

    def test_IO(self):
        pass
        #self.assertEquals(len(self.bigg), 323)

    def test_2(self):
        #i=0
        bigg_list=self.bigg
        #print(bigg_list[i])
        #bigg_list[i].to_pyrange(0)
        #print(i, bigg_list[i].ExonPyrange)
        #print(i, bigg_list[i].IntronPyrange)

        for j in range(2, 5):
            print(bigg_list[j])
            bigg_list[j].to_pyrange(gene_start=0)
            print bigg_list[j].exon
            #print(j,bigg_list[j].ExonPyrange)
            #print(j,bigg_list[j].IntronPyrange)

        gene_start=0
        #distance_exon = bigg_list[i].pyrange_cal_distance_exon(bigg_list[j], gene_start, by="ratio")

    def test_cluster(self):
        pass
        sample = self.bigg[0:40]
        D=cal_distance(sample, intronweight=0.5,by="ratio", core=40)
        print(D)
        #keep=filter_D(D, sample)
        #print len(bigg_list)

    def tearDown(self):
        self.bigg = None
