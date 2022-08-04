#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/27/2018 9:18 AM
# @Author  : Runsheng     
# @File    : tracklist_test.py


from trackcluster.tracklist import *
from trackcluster.post import *
import unittest
import logging


class TracklistTest(unittest.TestCase):
    def setUp(self):
        self.biggfile="./genes/unc52/unc52_sw.bed"
        self.swexonfile="./genes/unc52/unc52_sw_exon.bed"
        self.csvfile="./genes/unc52/exon_inter.bed"
        #self.swfile="/home/zhaolab1/data/nanorna/score_SL1_ssw.txt"
        self.testout="./genes/unc52/test.bed"
        self.test_interout="./genes/unc52/exon_inter.bed"

    def test_read(self):
        self.bigg_list=read_bigg(self.biggfile)
        print(len(self.bigg_list))

    def test_boundary(self):
        self.bigg_test=read_bigg(self.testout)
        bigg0=self.bigg_test[0]
        bigg1=self.bigg_test[-1]

        #matched=boundary_compare(bigg0, bigg1)
        #print matched,

    def test_boundaryall(self):
        self.bigg_test=read_bigg(self.testout)
        list_ref=[x for x in self.bigg_test if x.ttype!="nanopore_read"]

        for i in list_ref:
            i.get_junction()

        for bigg in self.bigg_test:
                is_new=has_new_junction(bigg, list_ref)
                #print(bigg.name, is_new)
        logging.info("Test {} passed.".format(self.test_boundaryall.__name__))

    def test_class4(self):
        self.bigg_test=read_bigg(self.testout)
        list_ref=[x for x in self.bigg_test if x.ttype!="nanopore_read"]

        names=["930b8415-579c-4f85-bdfa-c0ddbf49719b","09b8ff24-e1bd-40de-9b47-1412ced6500e","eb6e13ea-56cb-437c-8aa2-a5140590d4e6","7be0aa03-ff5c-4e81-8e3f-ea329f2d2ac0","9c148ee0-1108-43c9-a7b5-6545d1922a92","c74e83ed-2e4e-4b89-b2e7-ca107fe8d73b","4f9d09ff-7bf6-46ca-b275-75553967ca74","8ef6b155-5fa9-4d84-b129-6794a73756bd","61efb881-efc4-49a9-a07f-9ba0eceb9a09","1a5c37df-a572-4d4b-90d3-b51d8df73e26","765ea80b-7b8e-42bf-b153-fc6d7d64f197","39e224fe-394a-4b98-9759-c2271fdb8799","7545fe38-5d0b-4e89-ad47-929d08d70282","48afbbef-bbfb-4915-9bf9-5705339e491f"]
        for bigg in self.bigg_test:
            if bigg.name in names[0:]:
                class4=class_4(bigg, list_ref, 10)
                print(bigg.name, class4)

    def test_pandas_summary(self):
        csvfile=self.csvfile
        i_dic=pandas_summary(csvfile)
        print(("length",len(i_dic)))

    def test_IO(self):
        pass
        #add_sw(self.biggfile, self.swfile, out="./test/unc52_sw.bed")

    def test_add_subread_bigg(self):
        self.bigg_list=read_bigg(self.biggfile)
        bigg_add=self.bigg_list+self.bigg_list+self.bigg_list
        bigg_added=add_subread_bigg(bigg_add)

        print(len(self.bigg_list)==len(bigg_added))

        for i,j in zip(self.bigg_list, bigg_added):
            if i.name==j.name:
                #print("Pass")
                pass
            else:
                print("Not equal in", i, j)

    def test_tobedfile(self):
        dir= "./genes/unc52"
        prefix=self.biggfile.split("/")[-1].split(".")[0]
        print(prefix)

        self.bigg_list=read_bigg(self.biggfile)
        bigglist_to_bedfile(self.bigg_list, dir=dir, prefix=prefix)

    def test_wrapper_bedtools(self):
        bed1=self.swexonfile
        out=wrapper_bedtools_intersect2(bed1, bed1, self.test_interout)
        from trackcluster.utils import count_file
        print(("file line count",  out , count_file(out)))

    def test_bigglist_add(self):
        pass

    def test_cat_bed(self):
        import os
        wkdir="./genes"
        os.chdir(wkdir)
        cat_bed( "**/*_simple_coveragej.bed")


    def tearDown(self):
        self=None
