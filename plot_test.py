import unittest

from track import bigGenePred
from plots import *

class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, True)


class PlotTest(unittest.TestCase):
    def setUp(self):
        bigg_nano=[]
        with open("./test/genes/unc52/unc52_sw.bed") as f:
            for line_one in f.readlines():
                bigg_one=bigGenePred()
                bigg_one.from_string(line_one)
                bigg_nano.append(bigg_one)
        bigg_gff=[]
        with open("./test/genes/unc52/unc52_gff.bed") as f:
            for line_one in f.readlines():
                bigg_one=bigGenePred()
                bigg_one.from_string(line_one)
                bigg_gff.append(bigg_one)
        self.bigg_nano=bigg_nano
        self.bigg_gff=bigg_gff

    def test_IO(self):
        self.assertEquals(len(self.gff),17)

    def test_plot(self):
        sample=self.bigg_nano[0:]
        line_plot_merge(sample, self.bigg_gff,
                        out="./test/genes/unc52/bb.pdf",
                  biggout="./test/genes/unc52/test.bed",
                  Dout="./test/genes/unc52/d.csv",
                  intronweight=0.5,
                  by="ratio_all")

    def test_plot_tree(self):
        pass
        #sample=self.bigg[0:10]
        #plot_tree(sample, out="./test/aa.pdf")

    def tearDown(self):
        self.bigg_nano=None
        self.bigg_gff=None


if __name__ == '__main__':
    unittest.main()
