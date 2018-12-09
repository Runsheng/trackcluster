#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/27/2018 9:18 AM
# @Author  : Runsheng     
# @File    : tracklist_test.py


from tracklist import *
import unittest


class TracklistTest(unittest.TestCase):
    def setUp(self):
        self.biggfile="./test/unc52_sw.bed"
        self.swfile="/home/zhaolab1/data/nanorna/score_SL1_ssw.txt"


    def test_read(self):
        self.bigg_list=read_bigg(self.biggfile)
        print len(self.bigg_list)

    def test_IO(self):
        pass
        #add_sw(self.biggfile, self.swfile, out="./test/unc52_sw.bed")

    def test_add_subread_bigg(self):
        self.bigg_list=read_bigg(self.biggfile)
        bigg_add=self.bigg_list+self.bigg_list+self.bigg_list
        bigg_added=add_subread_bigg(bigg_add)

        print len(self.bigg_list)==len(bigg_added)

        for i,j in zip(self.bigg_list, bigg_added):
            if i.name==j.name:
                pass
            else:
                print "Not equal in", i, j



    def test_tobedfile(self):
        dir="./test"
        prefix=self.biggfile.split("/")[-1].split(".")[0]
        print prefix

        self.bigg_list=read_bigg(self.biggfile)
        bigglist_to_bedfile(self.bigg_list, dir=dir, prefix=prefix)

    def test_wrapper_bedtools(self):
        bed1="./test/unc52_sw_exon.bed"
        out=wrapper_bedtools_intersect2(bed1, bed1, "./test/exon_inter.bed")
        print out

    def test_bigglist_add(self):
        pass

    def tearDown(self):
        pass
