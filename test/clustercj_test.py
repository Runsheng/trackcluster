#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2/13/2021 1:46 PM
# @Author  : Runsheng     
# @File    : clustercj_test.py
"""
Test set for clustercj mode
"""

from trackcluster.clustercj import *
from trackcluster.tracklist import *
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

        ### used the junction_pre to filter the reads first
        self.bigg_nano=junction_pre(self.bigg_nano, self.bigg_gff)
        self.bigg=self.bigg_nano+bigg_gff

        self.bigg_dic=list_to_dic(self.bigg)


    def test_get_junction_dic(self):
        self.site_dic=get_junction_dic(self.bigg)
        print(self.site_dic)

    def test_generate_binary_junction_list(self):
        self.site_dic=get_junction_dic(self.bigg)
        keys=list(self.site_dic.keys())
        print(keys)
        junction=self.bigg[10].junction
        print((generate_binary_junction_list(junction, keys)))

    def __test_get_read_junction_D(self):
        # test if using pandas dataframe is better than using dict with array as values
        self.site_dic=get_junction_dic(self.bigg)
        df=get_read_junction_dic(self.bigg, self.site_dic)
        #print(df)
        print (df.index.names)
        print(df["087b20ed-78ba-48f9-a038-f23a9c4b75c6", 14647857]==1)
        #df = df[df.duplicated(keep=False)]
        #print df

        dup_l=df.groupby(df.columns.tolist()).apply(lambda x: tuple(x.index)).tolist()
        n=0
        for i in dup_l:
            n+=len(i)
        print(("{} tracks have identical junctions".format(n)))


    def test_get_read_junction_dic(self):
        self.site_dic=get_junction_dic(self.bigg)
        df=get_read_junction_dic(self.bigg, self.site_dic)
        #print(df.keys())
        # only work for unc-52
        #print numpy.array(df.loc["087b20ed-78ba-48f9-a038-f23a9c4b75c6"])-numpy.array(df.loc["ZC101.2g.1"])
        print((df["1f940a05-24a1-4455-a506-b1aa04caf81a"]))
        print((df["087b20ed-78ba-48f9-a038-f23a9c4b75c6"]))
        print((df["1f940a05-24a1-4455-a506-b1aa04caf81a"]-df["087b20ed-78ba-48f9-a038-f23a9c4b75c6"]))
        print((df["087b20ed-78ba-48f9-a038-f23a9c4b75c6"]-df["1f940a05-24a1-4455-a506-b1aa04caf81a"]))

    def test_get_arrry_freq_is_junction_inside(self):
        self.site_dic=get_junction_dic(self.bigg)
        df=get_read_junction_dic(self.bigg, self.site_dic)
        #j1=df["087b20ed-78ba-48f9-a038-f23a9c4b75c6"]
        #j2=df["1f940a05-24a1-4455-a506-b1aa04caf81a"]
        #print(array_del)
        #print(get_array_freq(array_del))
        ### test is _junction_inside
        #print(is_junction_inside(j1, j2))
        #print(is_junction_inside(j2, j1))


        keys=list(df.keys())

        # test different method for speed
        #####-----
        # take 4s for unc-52, for reducing the caculation to 50%
        def iter_for():
            import itertools
            for i,j in itertools.combinations(keys, 2):
                _is_junction_inside(df[i], df[j])
        iter_for()
        #####-----
        # takes 8s
        def mutual_for():
            for i in keys:
                for j in keys:
                    _is_junction_inside(df[i], df[j])
        #mutual_for()

        # map method, 4s for unc52, 19s for

        def map_two():
            import itertools
            ij_l=itertools.combinations(keys,2)
            out=list(map(lambda x: _is_junction_inside(df[x[0]], df[x[1]]),
                         ij_l))
            print(len(out))
        #map_two()


    def __test_compare_junction(self):
        self.site_dic=get_junction_dic(self.bigg)
        df=get_read_junction_dic(self.bigg, self.site_dic)
        j_all = []
        for i in self.bigg:
            i.get_junction()
            for j in self.bigg[2:]:
                j.get_junction()
                j_all.append(__compare_junction(i.junction, j.junction))
        print(j_all)

    def test_get_corrected_dic(self):
        self.site_dic=get_junction_dic(self.bigg)
        get_corrected_dic(self.site_dic, 2,10)
        get_corrected_dic(self.site_dic, 3,10)
        get_corrected_dic(self.site_dic, 5,10)

    def test_bigg_correct(self):
        import copy
        self.site_dic=get_junction_dic(self.bigg)
        w_to_r, w_to_no=get_corrected_dic(self.site_dic, 2,20)
        count=0
        for bigg in self.bigg:
            bigg.bk=copy.deepcopy(bigg)
            bigg_1, flag=get_bigg_correct(bigg, w_to_r)
            if flag==1:
                #print bigg_1
                #print bigg.bk
                count+=1
        print(( "The number of reads corrected is ", count))

    def test_flow_junction_correct(self):
        bigg_correct, bigg_rare= flow_junction_correct(self.bigg, 2, 10)
        print((len(bigg_correct), len(set(bigg_correct)), len(bigg_rare), len(set(bigg_rare)) ))

    def test_flow_clusterj_corrected(self):
        bigg_correct, bigg_rare= flow_junction_correct(self.bigg)




    def __test_group_nearby_site(self):
        sites=[
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 20, 21, 25, 26, 27, 28, 37, 38, 39, 40],
            [0, 1, 2, 3, 4, 5,  9, 10, 11, 12, 13, 17,  25, 26, 27, 28, 37, 38, 39, 40],
            [],
            [0,1,9],
            [0,1,8,9],
            [14647306, 14647307, 14647857, 14647919, 14648142, 14648551, 14648827]
        ]

        for site in sites:
            print(__group_nearby_site(site))

        # expected out
        #[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], [20, 21, 25, 26, 27, 28, 37, 38, 39, 40]]
        #[[0, 1, 2, 3, 4, 5, 9, 10, 11, 12, 13], [17, 25, 26, 27, 28, 37, 38, 39, 40]]
        #[]
        #[[0, 1, 3]]





