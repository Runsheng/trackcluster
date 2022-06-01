#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/9/2018 10:53 AM
# @Author  : Runsheng     
# @File    : cluster_test.py

from trackcluster.clusterj import *
from trackcluster.tracklist import *
from trackcluster.tracklist import write_bigg
import unittest


class ClusterjTest(unittest.TestCase):
    def setUp(self):
        genes=["unc52", "AT1G06860", "AT2G02100", "AT2G43410"]
        gene=genes[0]

        bigg_nano=[]
        with open("./genes/{gene}/{gene}_nano.bed".format(gene=gene)) as f:
            for line_one in f.readlines():
                bigg_one=bigGenePred()
                bigg_one.from_string(line_one)
                bigg_nano.append(bigg_one)

        bigg_gff=[]
        with open("./genes/{gene}/{gene}_gff.bed".format(gene=gene)) as f:
            for line_one in f.readlines():
                bigg_one=bigGenePred()
                bigg_one.from_string(line_one)
                bigg_gff.append(bigg_one)

        self.bigg_nano=bigg_nano
        self.bigg_gff=bigg_gff
        self.bigg=bigg_nano+bigg_gff


    def test_get_corrected_junction(self):
        import random
        from copy import deepcopy
        #bigg=self.bigg[0]
        bigg=chose_read_from_list("1d7aab69-44a6-461d-81e3-e0b9e1269bdc", self.bigg)
        bigg.get_junction()
        print(bigg)
        exon_bk=deepcopy(bigg.exon)

        # mask some wrong junctions for the bigg
        junction_mask=[]
        for i in bigg.junction:
            offset_l=[-5,-4, -3,-2,-1,0,1,2,3,4,5]
            offset=random.choice(offset_l)
            junction_mask.append(i+offset)

        bigg.junction=junction_mask
        bigg.write_junction_to_exon()
        bigg.exon_to_block()
        print(bigg)
        ####

        junction_dic=get_junction_dic(self.bigg)

        junction_new=get_corrected_junction(bigg, junction_dic,
                         coverage_cutoff=2, offset=5)

        bigg.junction=junction_new
        bigg.write_junction_to_exon()
        bigg.exon_to_block()
        print(bigg)
        print(bigg.exon==exon_bk)


    def test_is_bigg1_inside_bigg2_junction(self):
        bigg1, bigg2 = (self.bigg[0], self.bigg[1])
        print(bigg1)
        print(bigg2)

        print(is_junction_inside(bigg1, bigg2))
        print(is_junction_inside(bigg2, bigg1))

    def test_junction_simple_merge(self):

        print(len(self.bigg))
        bigg_n=junction_simple_merge(self.bigg)

        print(len(bigg_n))
        write_bigg(bigg_n, "./test/genes/AT2G43410/tt.bed")

    def test_is_single_exon_in(self):
        # get single and note single
        single=[]
        muti=[]

        for bigg in self.bigg:
            bigg.get_junction()
            if len(bigg.junction)==0:
                single.append(bigg)
            else:
                muti.append(bigg)

        for s1 in single:
            for s2 in muti:
                print(s1)
                print(s2)
                print(is_single_exon_in(s1, s2))

    def test_flow_junction_cluster(self):
        bigg_subread=flow_junction_cluster(self.bigg_nano, self.bigg_gff)
        for i in bigg_subread:
            print(i)

        #write_bigg(bigg_subread, "./test/genes/AT2G43410/AT2G43410_subread.bed")

    def test_bigg_get_namedic(self):
        bigg_list=read_bigg("genes/AT2G43410/AT2G43410_subread.bed")
        name_dic=bigg_get_namedic(bigg_list)
        print(name_dic)

