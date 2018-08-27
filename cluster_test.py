#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/9/2018 10:53 AM
# @Author  : Runsheng     
# @File    : cluster_test.py

from cluster import *
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

    def test_cluster(self):
        sample = self.bigg[30:80]
        D,bigg_list=cal_distance(sample, filter=True,intronweight=0.5,by="ratio", core=40)
        #print(D)
        #keep=filter_D(D, sample)
        print len(bigg_list)


    def tearDown(self):
        self.bigg = None
