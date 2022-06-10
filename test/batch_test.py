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

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
