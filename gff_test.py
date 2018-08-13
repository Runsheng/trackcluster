#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/10/2018 10:45 AM
# @Author  : Runsheng     
# @File    : gff_test.py


from gff import *
import unittest


class GFFTest(unittest.TestCase):
    def setUp(self):
        self.gff = GFF("./test/unc52.gff")

    def test_format_write(self):
        self.gff.gene_format()
        print(len(self.gff.gene_d.keys()))
        print(self.gff.gene_d.keys()[0:10])
        #self.gff.gff_write("./test/unc52_write.gff", keys=["ZC101.2"])

    def test_parser(self):
        print(parse_attributes("ID=Transcript:ZC101.2e;Parent=Gene:WBGene00006787;Name=ZC101.2e;wormpep=WP:CE18424;locus=unc-52"))

    def test_to_transcript(self):
        self.gff.transcript_format()
        print len(self.gff.transcript_d), len(self.gff.transcript_to_gene)
        #print(sorted(self.gff.transcript_to_gene.keys()))
        print(self.gff.transcript_d)
        print("Finish to_transcript test")


    def tearDown(self):
        self.bigg = None



