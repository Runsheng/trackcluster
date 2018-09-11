#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/1/2018 2:21 PM
# @Author  : Runsheng     
# @File    : utils_test.py

import unittest

from utils import *


class PlotTest(unittest.TestCase):
    def setUp(self):
        self.bedfile1="/run/user/1002/ZC101.2a.1_exon.bed"
        self.bedfile2="/run/user/1002/ZC101.2a.1_intron.bed"

        self.bedfile_list=["/run/user/1002/ZC101.2a.1_exon.bed",
                           "/run/user/1002/ZC101.2a.1_intron.bed"]

    def test_IO(self):
        out=wrapper_bedtools_jaccard(self.bedfile1, self.bedfile2)
        print out

    def test_inter(self):
        out=wrapper_bedtools_intersect(self.bedfile1, self.bedfile_list)
        print out




