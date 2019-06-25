#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/29/2018 2:34 PM
# @Author  : Runsheng     
# @File    : batch_test.py


"""
insert a full function key value test for the plotlib
"""

import unittest

from track import bigGenePred
from batch import *
from tracklist import *

import os

class BatchTest(unittest.TestCase):
    def setUp(self):
        pass

    def test_process_1(self):
        os.chdir("./test/genes/unc52")
        key="unc52"
        process_one_subsample(key, intronweight=0.5,batchsize=500, full=True)

    def test_process_2(self):
        os.chdir("./test/genes/")
        key="AT1G06860"
        process_one_subsample(key, intronweight=0.5,batchsize=1000, full=True)

    def test_re_cal(self):

        bigg_nano=read_bigg("./test/unc52/unc52_simple_coverage.bed")
        bigg_gff=[]
        flow_cluster(bigg_nano, bigg_gff)


    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
