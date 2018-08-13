import unittest

from track import bigGenePred
from plots import *

class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, True)


class PlotTest(unittest.TestCase):
    def setUp(self):
        bigg=[]
        with open("./test/unc52.bed") as f:
            for line_one in f.readlines():
                bigg_one=bigGenePred()
                bigg_one.from_string(line_one)
                bigg.append(bigg_one)
        with open("./test/unc52_gff.bed") as f:
            for line_one in f.readlines():
                bigg_one=bigGenePred()
                bigg_one.from_string(line_one)
                bigg.append(bigg_one)
        self.bigg=bigg

    def test_IO(self):
        self.assertEquals(len(self.bigg),339)

    def test_plot(self):
        sample=self.bigg[0:]
        line_plot(sample, out="./test/bb.pdf")

    def test_plot_tree(self):
        pass
        #sample=self.bigg[0:10]
        #plot_tree(sample, out="./test/aa.pdf")

    def tearDown(self):
        self.bigg=None


if __name__ == '__main__':
    unittest.main()
