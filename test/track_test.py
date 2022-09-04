#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/1/2018 12:03 AM
# @Author  : Runsheng     
# @File    : track_test.py


import unittest

from trackcluster.track import *
from trackcluster.utils import fasta2dic

class Test(unittest.TestCase):
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

    def test_IO(self):
        #sample=self.bigg[0]
       for sample in self.bigg:
            sample.to_bedstr()
            # test single exon track
            if sample.intronlen==0:
                print(sample)
                print((sample.exon_str))
                print("======")
                print(sample.intron_str)

    def test_get_exon(self):
        sample = self.bigg[0]
        if sample.name=="579ebc2e-86ca-469f-b68c-7262fc292c9d":
            sample.get_exon()

            print((sample.exon==[(14627636, 14627786), (14627936, 14628112), (14628157, 14628269), (14628325, 14628427), (14628724, 14630134), (14649328, 14649363)]))
            print((sample.intron==[(14627786, 14627936), (14628112, 14628157), (14628269, 14628325), (14628427, 14628724), (14630134, 14649328)]))
            print((sample.exonlen==1985, sample.intronlen==19742))
        else:
            print("Not run test_get_exon")

    def test_write_junction_to_exon(self):
        from copy import deepcopy
        ### need to test one foward and one reverse isoform
        # use unc-52 and AT2G43410 in init

        sample = self.bigg[0]
        sample.get_junction()
        print((sample.junction))
        print((sample.exon))
        sample_bk=deepcopy(sample.exon)
        # re-init
        sample.exon=None
        sample.intron=None
        sample.exonlen=0
        sample.intronlen=0

        sample.write_junction_to_exon()
        print((sample.exon))
        print((sample.exon==sample_bk))

        # for unc-52 only
        #print sample.exonlen==1985, sample.intronlen==19742

    def test_exon_to_block(self):
        sample = self.bigg[0]
        print(sample)
        print((sample.chromStarts, sample.blockSizes))
        sample.get_junction()

        # re-init
        sample.chromStarts=[]
        sample.blockSizes=[]
        print(sample)

        sample.exon_to_block()
        print(sample)

    def test_cal_distance(self):

        bed1, bed2=self.bigg[0:2]
        #print(bed1.bedfile_cal_distance_exon(bed2))
        #print(bed1.bedfile_cal_distance_intron(bed2))
        pass


    def __test_bindseq_reverse(self):
        ### need to use one gene
        ref_dict=fasta2dic("/home/li/reference/tair/tair10.fa")
        mrna_dict=fasta2dic("./test/genes/AT2G43410/AT2G43410.1.fasta")

        name="rna-NM_001337028.1"
        bigg_one=None
        for i in self.bigg:
            if i.name==name:
                bigg_one=i
        if bigg_one is None:
            return False
        print((bigg_one.name))

        real_seq=str(mrna_dict[name].seq)

        bigg_one.bind_chroseq(ref_dict, gap=0, intron=False)
        real_seq=real_seq.replace("\n", "")
        print((bigg_one.seq_chro==real_seq))
        print(bigg_one)

        pos=bigg_one.mrna_pos_to_chro(0)
        print(pos)

        # for single site
        print((chr_select(ref_dict, "chr2", 18026397, 18026398)))
        print((chr_select(ref_dict, pos[0], pos[1], pos[1]+1 )))

    def test_bindseq_forward(self):
        #ref_dict = fasta2dic("/home/li/reference/tair/tair10.fa")
        mrna_dict = fasta2dic("./genes/AT2G02100/AT2G02100.fasta")


    def test_mrna_pos_to_chro(self):
        sample = self.bigg[0]
        print(sample)
        print((sample.mrna_pos_to_chro(100),
        sample.mrna_pos_to_chro(1000),
        sample.mrna_pos_to_chro(-1)))

        ### check if the nucl of the same pos is indentical
        ### chro is 0 based and mRNA is 1 based

    def test_orfs(self):
        """
        can only run with ce10 ref
        """
        # ref_dic is a large external file
        ref_dict=fasta2dic("/data/reference/ce10_ucsc.fa")
        bigg_one=self.bigg[30]
        bigg_one.bind_chroseq(ref_dict, gap=0, intron=False)
        print((bigg_one.seq_chro))
        ans=bigg_one.find_orfs_with_trans()

        print(ans)
        print(bigg_one)

    def test_orf(self):
        pass


if __name__ == '__main__':
    unittest.main()




