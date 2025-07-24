#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 4/11/2025 12:13 AM
# @Author  : Runsheng     
# @File    : fusion_test.py



from trackcluster.tracklist import read_bigg, write_bigg
from trackcluster.fusion import *
import unittest

from trackcluster.utils import group_site


class TracklistTest(unittest.TestCase):
    def setUp(self):
        self.prefix="cel"
        self.wkdir="/t1/celtrack"
        self.ref="refs.bed"
        self.read="reads_sorted.bed"

    def test_flow_fusion(self):
        flow_fusion(wkdir=self.wkdir,
                    prefix=self.prefix,
                    bigg_gff_file=self.ref,
                    bigg_nano_file=self.read,
                    f1=0.1, f2=0.1)

        flow_fusion(wkdir=self.wkdir,
                    prefix=self.prefix,
                    bigg_gff_file=self.ref,
                    bigg_nano_file=self.read,
                    f1=0.2, f2=0.2)

        flow_fusion(wkdir=self.wkdir,
                    prefix=self.prefix,
                    bigg_gff_file=self.ref,
                    bigg_nano_file=self.read,
                    f1=0.3, f2=0.3)











