from trackcluster.tracklist import *
from trackcluster.pre import *
import unittest
import logging


class PreTest(unittest.TestCase):
    def setUp(self):
        genes=["unc52", "AT1G06860", "AT2G02100", "AT2G43410"]
        gene=genes[3]

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

        self.bigg_nano_name="./genes/{gene}/{gene}_nano.bed".format(gene=gene)
        self.bigg_gff_name="./genes/{gene}/{gene}_gff.bed".format(gene=gene)


        self.merge_out_name="./genes/{gene}/{gene}_merge.bed".format(gene=gene)


    def test_wrapper_bedtools_intersect(self):
        outfile=wrapper_bedtools_intersect2_select(self.bigg_nano_name, self.bigg_nano_name)
        print(outfile)

    def test_group_bigg_by_gene(self):

        gene_bigg, bigg_novel_gene= group_bigg_by_gene(self.bigg_gff)
        print(gene_bigg)
        print(bigg_novel_gene==[]) # should be empty

    def test_wrapper_bedtools_merge(self):
        out=wrapper_bedtools_merge(self.bigg_nano_name, out=self.merge_out_name)
        print(out)











    def tearDown(self) -> None:
        self=None
