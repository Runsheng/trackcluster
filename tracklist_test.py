#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/27/2018 9:18 AM
# @Author  : Runsheng     
# @File    : tracklist_test.py


from tracklist import *
import unittest


class TracklistTest(unittest.TestCase):
    def setUp(self):
        self.biggfile="./test/unc52.bed"
        self.swfile="/home/zhaolab1/data/nanorna/score_SL1_ssw.txt"


    def test_IO(self):

        add_sw(self.biggfile, self.swfile, out="./test/unc52_sw.bed")


    def tearDown(self):
        pass
