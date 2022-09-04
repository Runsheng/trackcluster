#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 13/11/2018 6:52 PM
# @Author  : Runsheng
# @File    : utils_test.py.py


import unittest
from trackcluster.utils import *

class myUtilsTest(unittest.TestCase):

    def setUp(self):
        self.fafile="./genes/unc52/test_io.fa"

    def test_chro_select(self):
        seq_dic=fasta2dic(self.fafile)
        name, seq=chr_select(seq_dic, "test", 0, 12)
        print(seq)

        seq_rc=reverse_complement(seq)
        print(seq_rc)

    def test_myexe(self):
        #print(myexe("ls -a |grep init"))
        print((myexe("bedtools --version")))
        print(myexe("which bedtools"))
        print("minimap2 --version")
        print(myexe("echo $PATH"))

    def test_is_bin_in(self):
        print(is_bin_in("minimap2"))
        print(myexe("minimap2 --help"))
        import shutil
        print(shutil.which("minimap2"))

        myexe("/bin/bash -i -c minimap2")

    def test_is_package_installed(self):
        print(is_package_installed("Bio"))

    def test_summary(self):
        logger = log_summary()
        logger.info("This is just a test, print number 100 as float" + "," + str(float(1) * 100))
        logger.info("In total %d line, %s" % (2, "just a test number"))

    def test_detail_file(self):
        logger = log_detail_file("./log.txt")
        logger.debug("This is just a test, print number 100 as float" + "," + str(float(1) * 100))
        logger.debug("In total %d line, %s" % (2, "just a test number"))

    def tearDown(self):
        self.fafile=None

if __name__ == '__main__':
    unittest.main()