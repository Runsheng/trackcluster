#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2/20/2019 12:13 AM
# @Author  : Runsheng     
# @File    : post_test.py


#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/27/2018 9:18 AM
# @Author  : Runsheng
# @File    : tracklist_test.py

from trackcluster.tracklist import read_bigg, write_bigg
from trackcluster.post import *
import unittest

from trackcluster.utils import group_site


class TracklistTest(unittest.TestCase):
    def setUp(self):
        self.biggfile="./genes/unc52/unc52_sw.bed"
        #self.swfile="/home/zhaolab1/data/nanorna/score_SL1_ssw.txt"
        self.testout="./genes/unc52/test.bed"
        self.name=["930b8415-579c-4f85-bdfa-c0ddbf49719b","09b8ff24-e1bd-40de-9b47-1412ced6500e","eb6e13ea-56cb-437c-8aa2-a5140590d4e6","7be0aa03-ff5c-4e81-8e3f-ea329f2d2ac0","9c148ee0-1108-43c9-a7b5-6545d1922a92","c74e83ed-2e4e-4b89-b2e7-ca107fe8d73b","4f9d09ff-7bf6-46ca-b275-75553967ca74","8ef6b155-5fa9-4d84-b129-6794a73756bd","61efb881-efc4-49a9-a07f-9ba0eceb9a09","1a5c37df-a572-4d4b-90d3-b51d8df73e26","765ea80b-7b8e-42bf-b153-fc6d7d64f197","39e224fe-394a-4b98-9759-c2271fdb8799","7545fe38-5d0b-4e89-ad47-929d08d70282","48afbbef-bbfb-4915-9bf9-5705339e491f"]
        self.bigg_test=read_bigg(self.testout)

    def test_correction(self):
        self.bigg_test=read_bigg(self.testout)
        list_ref=[x for x in self.bigg_test if x.ttype!="nanopore_read"]

        pass


    def test_boundaryall(self):
        self.bigg_test=read_bigg(self.testout)
        list_ref=[x for x in self.bigg_test if x.ttype!="nanopore_read"]

        for i in list_ref:
            i.get_junction()

        for bigg in self.bigg_test:
                is_new=has_new_junction(bigg, list_ref)
                print(bigg.name, is_new)

    def test_class4(self):
        self.bigg_test=read_bigg(self.testout)
        list_ref=[x for x in self.bigg_test if x.ttype!="nanopore_read"]

        names=self.name
        for bigg in self.bigg_test:
            if bigg.name in names[0:]:
                class4=class_4(bigg, list_ref, 10)
                print(bigg.name, class4)

    def test_compare_ei_by_boudary(self):
        self.bigg_test=read_bigg(self.testout)
        list_ref=[x for x in self.bigg_test if x.ttype!="nanopore_read"]

        bigg5=[]
        for bigg in self.bigg_test:
            if bigg.name in self.name:
                bigg5.append(bigg)

        for bigg0 in bigg5:
            for bigg_ref in list_ref:
                m,e=compare_ei_by_boudary(bigg0, bigg_ref)
                if len(m)%2!=0 or len(e)%2!=0:
                    print(bigg0.name, bigg_ref.name, m, len(m), e, len(e))

    def test_des_ei_by_boudary(self):
        """
        use ZC101.2f.1 and  930b8415-579c-4f85-bdfa-c0ddbf49719b
        """
        list_used=[x for x in self.bigg_test if x.name in ["ZC101.2f.1",  "930b8415-579c-4f85-bdfa-c0ddbf49719b"]]
        list_used=[x for x in self.bigg_test if x.name in ["ZC101.2r.1",  "48afbbef-bbfb-4915-9bf9-5705339e491f",
                                                           "7545fe38-5d0b-4e89-ad47-929d08d70282"]]

        desc_ei_by_boundary(list_used[0], list_used[2])

    def test_find_nearest(self):

        self.bigg_test=read_bigg(self.testout)
        list_ref=[x for x in self.bigg_test if x.ttype!="nanopore_read"]

        for bigg in self.bigg_test:
            if bigg.name ==self.name[0]:
                bigg0=bigg

        nearest_ref=find_nearest_ref(bigg0, list_ref)
        bigg_l=[bigg0, nearest_ref]

        # write the bigg to check them manualy
        write_bigg(bigg_l, "./bigg2.bed")

    def test_flow(self):

        list_ref=[x for x in self.bigg_test if x.ttype!="nanopore_read"]

        bigg5=[]
        for bigg in self.bigg_test:
            if bigg.name in self.name:
                bigg5.append(bigg)

        print(len(bigg5))
        for bigg in bigg5:

            print(flow_desc(bigg, list_ref, offset=20))

    def test_group_site(self):

        missed=[
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 20, 21, 25, 26, 27, 28, 37, 38, 39, 40],
            [],
            [0,1,3]]

        for missed_one in missed:
            print(group_site(missed_one))


    def test_desc_group(self):
        groups= [[0, 1, 2, 3, 4, 5, 6, 7], [14, 15], [27, 28, 29, 30]]



    def test_validall(self):
        pass









