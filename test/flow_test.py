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
        flow_key_clusterj(wkdir, genename_file, 30)

    def test_mergebed(self):
        wkdir="/t1/shoudong_488/test/trackall"
        prefix="488_aba_1"
        novel_bed="488_aba_1_novel.bed"
        outfile=prefix+"_merge.bed"

        os.chdir(wkdir)
        out=wrapper_bedtools_merge(novel_bed, outfile)
        print(out)

    def test_mergebed2bigg(self):
        wkdir="/t1/shoudong_488/test/trackall"
        prefix="488_aba_1"
        novel_bed="488_aba_1_novel.bed"
        mergefile=prefix+"_merge.bed"
        bigg_regionmark=prefix+"_regionmark.bed"

        os.chdir(wkdir)
        bigg_l=mergedbed2bigg(mergefile, count_cutoff=5)

        write_bigg(bigg_l, bigg_regionmark)


    def test_wrapper_bedtools_subtract(self):

        wkdir="/t1/shoudong_488/test/trackall"
        prefix="488_aba_1"
        novel_bed="488_aba_1_novel.bed"
        bigg_regionmark="488_aba_1_regionmark.bed"
        gff_bed="gene.bed_s"
        bigg_regionmark_f=prefix+"_regionmarkf.bed"

        os.chdir(wkdir)

        # 421 for the new genes,
        bigg_regionmark_f=wrapper_bedtools_subtract(bigg_regionmark, gff_bed,bigg_regionmark_f)


    def test_prepar_noveldir(self):

        wkdir="/t1/shoudong_488/test/trackall"
        prefix="488_aba_1"
        novel_bed="488_aba_1_novel.bed"
        bigg_regionmark="488_aba_1_regionmarkf.bed"
        novelname_file="novelname.txt"


        # need to re-run flow_preparedir test to make the result correct
        # the substract will change the original novel file
        flow_preparedir(wkdir,prefix, bigg_regionmark, novel_bed, genename_file=novelname_file)

        # run
        flow_key_clusterj(wkdir, novelname_file, core=30)

    def test_prepare_run_gene_novel(self):
        """
        batch test for full gene and novel run
        :return:
        """
        # parameters
        wkdir="/t1/shoudong_488/test/tracktest"
        prefix="488_aba_1"
        gff_bed="../gene.bed_s"
        nano_bed="../488_aba_1_s.bed"
        f1=0.01
        f2=0.05
        core=30
        os.chdir(wkdir)

        flow_clusterj_all_gene_novel(nano_bed=nano_bed, gff_bed=gff_bed, prefix=prefix,wkdir=wkdir,
                                     core=core, f1=f1, f2=f2)


    def test_prepare_run_gene_novel_cluster_original(self):
        """
        batch test for full gene and novel run
        :return:
        """
        # parameters
        wkdir="/t1/shoudong_488/test/tracktest_noj"
        prefix="488_aba_1"
        gff_bed="../gene.bed_s"
        nano_bed="../488_aba_1_s.bed"
        f1=0.01
        f2=0.05
        core=30
        batchsize=2000
        intronweight=0.2
        cutoff1=0.05
        cutoff2=0.05
        os.chdir(wkdir)

        flow_cluster_all_gene_novel(nano_bed=nano_bed, gff_bed=gff_bed, prefix=prefix,wkdir=wkdir,
                                     core=core, f1=f1, f2=f2, intronweight=intronweight, batchsize=batchsize,
                                    cutoff1=cutoff1, cutoff2=cutoff2)


    def test_count(self):
        wkdir = "/t1/shoudong_488/test/tracktest"
        prefix = "488_aba_1"
        isoform_bed="488_aba_1_isoform.bed"
        gff_bed = "../gene.bed_s"
        nano_bed = "../488_aba_1_s.bed"
        os.chdir(wkdir)
        flow_count(wkdir, prefix, nano_bed, isoform_bed, gff_bed)


    def tearDown(self) -> None:
        self=None

