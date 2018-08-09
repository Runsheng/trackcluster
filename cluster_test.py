#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/9/2018 10:53 AM
# @Author  : Runsheng     
# @File    : cluster_test.py

from cluster import *
import unittest


class ClusterTest(unittest.TestCase):
    def setUp(self):
        bigg = []
        with open("./test/unc52.bed") as f:
            for line_one in f.readlines():
                bigg_one = bigGenePred()
                bigg_one.from_string(line_one)
                bigg.append(bigg_one)
        self.bigg = bigg

    def test_IO(self):
        self.assertEquals(len(self.bigg), 323)

    def test_cluster(self):
        sample = self.bigg[0:50]
        D=cal_distance(sample)
        plot_cluster(D)


    def tearDown(self):
        self.bigg = None
