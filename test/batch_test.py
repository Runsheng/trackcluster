#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/29/2018 2:34 PM
# @Author  : Runsheng     
# @File    : batch_test.py


"""
insert a full function key value test for the plotlib
"""

import unittest

from trackcluster.batch import *

import os

from trackcluster.tracklist import cat_bed


class BatchTest(unittest.TestCase):
    def setUp(self):
        pass

    def test_process_1(self):
        os.chdir("genes")
        key="unc52"
        process_one_subsample_try(key, intronweight=0.5,batchsize=300, full=True)

    def test_process_2(self):
        os.chdir("genes")
        key="AT2G43410"
        process_one_subsample_try(key, intronweight=0.5,batchsize=1000, full=True)

    def test_process_3(self):
        # large reads number
        os.chdir("genes")
        key = "AT2G02100"
        #process_one_subsample_try(key, intronweight=0.5, batchsize=900, full=True)
        process_one_subsample_try(key, intronweight=0.5, batchsize=1000, full=True)  # 33S
        #process_one_subsample_try(key, intronweight=0.5, batchsize=2000, full=True) # 46S
        # should be 898

    def test_junction_3(self):
        os.chdir("genes")
        key = "AT2G02100"
        key = "unc52" # 90,40
        # size=2000,28s， size=1000, 13s; size=500,6s
        #out=process_one_junction_corrected_try(key, batchsize=2000, full=True)  # 27S, still slow
        #out=process_one_junction_corrected_try(key, batchsize=1000, full=True)  # slow

        out=process_one_junction_corrected_try(key, batchsize=200, full=True)  #  5s

    def test_cat_bed(self):
        wkdir="./genes"
        os.chdir(wkdir)
        cat_bed( "**/*_simple_coveragej.bed")

        wkdir="/t1/shoudong_488/test/tracktest"
        os.chdir(wkdir)
        cat_bed( "**/*_simple_coveragej.bed")



    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
