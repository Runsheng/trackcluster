import unittest

from track import bigGenePred
from plotsi import *

class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, True)


class PlotTest(unittest.TestCase):
    def setUp(self):
        bigg_nano=[]
        with open("./test/unc52_sw.bed") as f:
            for line_one in f.readlines():
                bigg_one=bigGenePred()
                bigg_one.from_string(line_one)
                bigg_nano.append(bigg_one)
        bigg_gff=[]
        with open("./test/unc52_gff.bed") as f:
            for line_one in f.readlines():
                bigg_one=bigGenePred()
                bigg_one.from_string(line_one)
                bigg_gff.append(bigg_one)
        self.bigg_nano=bigg_nano
        self.bigg_gff=bigg_gff

    def test_IO(self):
        self.assertEquals(len(self.bigg),340)

    def test_plot(self):
        sample=self.bigg_nano[0:10]
        line_plot_merge(sample, self.bigg_gff,
                        out="./test/bb.pdf",
                  biggout="./test/test.bed",
                  Dout="./test/d.csv",
                  intronweight=0.5,
                  by="ratio_all", core=40)

    def test_plot_tree(self):
        pass
        #sample=self.bigg[0:10]
        #plot_tree(sample, out="./test/aa.pdf")

    def tearDown(self):
        self.bigg_nano=None
        self.bigg_gff=None


if __name__ == '__main__':
    unittest.main()
