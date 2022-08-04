#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/10/2018 10:45 AM
# @Author  : Runsheng     
# @File    : gff_test.py


from trackcluster.gff import *
import unittest


class GFFTest(unittest.TestCase):
    def setUp(self):
        #self.gff = GFF("/home/zhaolab1/reference/ce10_ws266.gff")
        self.gff=GFF("./genes/unc52/unc52.gff")

    def test_format_write(self):
        self.gff.gene_format()
        print((len(list(self.gff.gene_d.keys()))))
        print((list(self.gff.gene_d.keys())[0:10]))
        self.gff.gff_write("./genes/unc52/unc52_write.gff", keys=["ZC101.2"])
        print("test gff IO done")

    def test_parser(self):
        print((parse_attributes("ID=Transcript:ZC101.2e;Parent=Gene:WBGene00006787;Name=ZC101.2e;wormpep=WP:CE18424;locus=unc-52")))

    def test_to_transcript(self):
        self.gff.transcript_format()
        print(len(self.gff.transcript_d), len(self.gff.transcript_to_gene))
        print((sorted(self.gff.transcript_to_gene.keys())))
        print((self.gff.transcript_d))
        print("Finish to_transcript test")


    def tearDown(self):
        self.bigg = None



