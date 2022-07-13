import unittest
from trackcluster.flow import *

class flowtest(unittest.TestCase):
    """
    The test case for flow is not included in test folder
    too large files for fastq and bam
    """
    def setUp(self):
        pass

    def test_flow_mappping(self):
        wkdir="/t1/shoudong_488/test"
        ref_file="c24.fasta"
        fastq_file="488_aba_1.fastq"
        prefix = get_file_prefix(fastq_file, sep=".")
        self.bam_name=flow_mapping(wkdir,ref_file, fastq_file,prefix=prefix,core=24)

    def test_flow_bamconvert(self):
        wkdir="/t1/shoudong_488/test"
        prefix="488_aba_1"
        bamfile="488_aba_1_s.bam"
        out=prefix+".bed"
        self.biggnano=flow_bamconvert(wkdir, bamfile, out=out, prefix=prefix, score=30)

    def test_prepare_dir(self):
        wkdir="/t1/shoudong_488/test/trackall"
        prefix="488_aba_1"
        gff_bed="gene.bed_s"
        nano_bed="488_aba_1_s.bed"

        # the novel bigg file is prefix_novel.bed
        # the gene name file for further use is prefix_gene.txt
        gene_l=flow_preparedir(wkdir, prefix, bigg_gff_file=gff_bed, bigg_nano_file=nano_bed)
        print(gene_l[100])


    def test_flow_gene(self):
        """
        use
        test the 488 test data
        :return:
        """
        wkdir="/t1/shoudong_488/test/trackall"
        prefix="488_aba_1"
        genename_file = prefix + "_gene.txt"

        # recover the gene name
        os.chdir(wkdir)
        flow_gene_clusterj(wkdir, genename_file, 30)



    def test_prepar_noveldir(self):
        wkdir="/t1/shoudong_488/test/trackall"
        prefix="488_aba_1"
        novel_bed="488_aba_1_novel.bed"















    def tearDown(self) -> None:
        self=None

