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
        self.attr_str="ID=Transcript:ZC101.2e;Parent=Gene:WBGene00006787;Name=ZC101.2e;wormpep=WP:CE18424;locus=unc-52"

    def test_format_write(self):
        self.gff.gene_format()
        print((len(list(self.gff.gene_d.keys()))))
        print((list(self.gff.gene_d.keys())[0:10]))
        self.gff.gff_write("./genes/unc52/unc52_write_2.gff")
        print("test gff IO done")

    def test_parser_attr(self):
        print(parse_attributes(self.attr_str))

    def test_attr2str(self):
        attr_dic=parse_attributes(self.attr_str)
        attr_str_new=attr_dic_to_str(attr_dic)
        print(attr_str_new)
        print(self.attr_str)

    def test_parser_gff_and_reverse(self):
        self.gff.to_list()
        gff_one=self.gff.gff_string_list[1]
        gff_record=parse_gff_line(gff_one)
        print(gff_record)
        gff_str=gff_record_to_str(gff_record)
        print(gff_one)
        print(gff_str)


    def test_to_transcript(self):
        self.gff.transcript_format()
        print(len(self.gff.transcript_d), len(self.gff.transcript_to_gene))
        print((sorted(self.gff.transcript_to_gene.keys())))
        print((self.gff.transcript_d))
        print("Finish to_transcript test")

    def test_ensembl(self):
        """
        test the gff performace in ensembl gff format
        :return:
        """
        ensembl_file="/data/reference/mouse/Mus_musculus.GRCm39.107.chr.gff3"
        gff=GFF(ensembl_file)
        gff.to_list()
        gff.gene_format()
        gff.transcript_format()
        print(len(gff.gene_d), len(gff.transcript_d))

        def __test_gene():
            """
            test code for ensembl mouse gff
            :return:
            """
            keys=["gene:ENSMUSG00000114538","gene:ENSMUSG00000103506" ]
            for k in keys:
                print(k)
                print(gene2line[k])
                for i in gene2line[k]:
                    print(self.gff_records[i])



    def test_ucsc(self):
        """
        todo: this need a gtg parser
        test the gff convert from ucsc gtf or gff
        :return:
        """
        pass


    def test_wormbase(self):
        pass

    def test_genbank(self):
        """
        test the gff format from genbank bioproject
        :return:
        """
        pass





    def tearDown(self):
        self.bigg = None



