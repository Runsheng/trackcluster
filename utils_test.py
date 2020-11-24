#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 13/11/2018 6:52 PM
# @Author  : Runsheng
# @File    : utils_test.py.py


import unittest
from utils import *

class myUtilsTest(unittest.TestCase):

    def setUp(self):
        self.fafile="./test/genes/unc52/test_io.fa"

    def test_chro_select(self):
        seq_dic=fasta2dic(self.fafile)
        name, seq=chr_select(seq_dic, "test", 0, 12)
        print(seq)

        seq_rc=reverse_complement(seq)
        print(seq_rc)


    def tearDown(self):
        self.fafile=None

if __name__ == '__main__':
    unittest.main()