#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/8/2018 2:44 PM
# @Author  : Runsheng     
# @File    : track.py

# third part import
from utils import set_tmp, wrapper_bedtools_jaccard
import os

dirp=set_tmp()

def get_start_end(list2):
    start = [x[0] for x in list2]
    end = [x[1] for x in list2]

    return start, end


class bigGenePred(object):
    """
    A data class for bigGenePred
    "bigGenePred gene models"
    (
    string chrom; "Reference sequence chromosome or scaffold"
    uint chromStart; "Start position in chromosome"
    uint chromEnd; "End position in chromosome"
    string name; "Name or ID of item, ideally both human readable and unique"
    uint score; "Score (0-1000)"
    char[1] strand; "+ or - for strand"
    uint thickStart; "Start of where display should be thick (start codon)"
    uint thickEnd; "End of where display should be thick (stop codon)"
    uint reserved; "RGB value (use R,G,B string in input file)"
    int blockCount; "Number of blocks"
    int[blockCount] blockSizes; "Comma separated list of block sizes"
    int[blockCount] chromStarts; "Start positions relative to chromStart"
    string name2; "Alternative/human readable name"
    string cdsStartStat; "Status of CDS start annotation (none, unknown, incomplete, or complete)"
    string cdsEndStat; "Status of CDS end annotation (none, unknown, incomplete, or complete)"
    int[blockCount] exonFrames; "Exon frame {0,1,2}, or -1 if no frame for exon"
    string type; "Transcript type"
    string geneName; "Primary identifier for gene"
    string geneName2; "Alternative/human readable gene name, used to store sample information such as stage or group"
    string geneType; "Gene type"
    )

    """
    def __init__(self):
        self.chrom=""
        self.chromStart=0
        self.chromEnd=0
        self.name=""
        self.score= 0
        self.strand="+"
        self.thickStart=0
        self.thickEnd=0
        self.reserved=[255,128,0]
        self.blockCount=0
        self.blockSizes=[]
        self.chromStarts=[]
        self.name2=""
        self.cdsStartStat="none"
        self.cdsEndStat="none"
        self.exonFrames=[-1]
        self.ttype="nanopore_read"
        self.geneName="none"
        self.geneName2="none"
        self.geneType="none"

        ### internal para
        self.exon=None
        self.intron=None
        self.interval_set=None
        self.gene_start=None
        self.exon_str=None
        self.intron_str=None
        self.exonlen=0
        self.intronlen=0

        self.exon_file=None
        self.intron_file=None

        self.seq=None
        self.subread=[] # use to store the reads contained inside the

    def to_list(self):
        data_l=[self.chrom,
                self.chromStart,
                self.chromEnd,
                self.name,
                self.score,
                self.strand,
                self.thickStart,
                self.thickEnd,
                ",".join([str(x) for x in self.reserved]),
                self.blockCount,
                ",".join([str(x) for x in self.blockSizes])+",",
                ",".join([str(x) for x in self.chromStarts])+",",
                self.name2,
                self.cdsStartStat,
                self.cdsEndStat,
                ",".join([str(x) for x in self.exonFrames])+",",
                self.ttype,
                self.geneName,
                self.geneName2,
                self.geneType]
        return data_l

    def to_str(self):
        """
        IO
        get a str in bigpred file
        """
        data_l=self.to_list()
        str_l=[str(x) for x in data_l]

        # ucsc compatiable chr name
        return "chr"+"\t".join(str_l)

    def __str__(self):
        return self.to_str()

    def from_string(self, string_bed):
        """
        IO
        reverse function of the to_str
        """
        string_l=string_bed.strip().split("\t")

        self.chrom=string_l[0][3:] # revsere the +chr operation
        self.chromStart=int(string_l[1])
        self.chromEnd=int(string_l[2])
        self.name=string_l[3]
        self.score=int(string_l[4])
        self.strand=string_l[5]
        self.thickStart=int(string_l[6])
        self.thickEnd=int(string_l[7])
        self.reserved=[int(x) for x in string_l[8].split(",")]
        self.blockCount=int(string_l[9])
        self.blockSizes=[int(x) for x in string_l[10][:-1].split(",")]
        self.chromStarts=[int(x) for x in string_l[11][:-1].split(",")]
        self.name2=string_l[12]
        self.cdsStartStat=string_l[13]
        self.cdsEndStat=string_l[14]
        self.exonFrames=[int(x) for x in string_l[15][:-1].split(",")]
        self.ttype=string_l[16]
        self.geneName=string_l[17]
        self.geneName2=string_l[18]
        self.geneType=string_l[19]

    def to_interval_set(self, gene_start=None):
        """
        for pair wise matrix generation
        :param gene_start:
        :return:
        """
        gene_start=self.chromStart if gene_start is None else gene_start
        if self.exon is None:
            self.get_exon(gene_start)

        interval_set=set()
        for i in self.exon:
            x,y=i
            set_one=set(range(x, x+y))
            interval_set=interval_set.union(set_one)
        self.interval_set=interval_set

    def __cal_distance(self, other_bgp, gene_start):
        """
        for pairwise compare
        compare the distance of this bgp with other bgp in same gene
        using the interval set
        """

        if self.interval_set is None:
            self.to_interval_set(gene_start)
        if other_bgp.interval_set is None:
            other_bgp.to_interval_set(gene_start)

        distance=len(self.interval_set.union(other_bgp.interval_set))-\
                 len(self.interval_set.intersection(other_bgp.interval_set))

        return distance

    def get_exon(self, gene_start=None):
        gene_start=self.chromStart if gene_start is None else gene_start
        offset=self.chromStart-gene_start

        line_exon=[]
        line_intron=[]

        len_count=0

        for x, y in zip(self.chromStarts, self.blockSizes):
            line_one=(x+offset, x+offset+y)
            line_exon.append(line_one)
            len_count+=y

        self.exon=line_exon
        self.exonlen=len_count

        for n, pair in enumerate(line_exon):
            start, end=pair
            if n==0:
                pass
            else:
                line_intron.append((end_p+1, start-1))
            if n==len(line_exon):
                break
            start_p, end_p=pair

        self.intron=line_intron
        len_intron=0
        for x,y in self.intron:
            len_intron+=(y-x)
        self.intronlen=len_intron

        # debug:
        #print("exon", self.exon)

    def to_bedstr(self, gene_start=None):
        """
        convert the exon region to BedTool object
        :return:
        """
        if self.exon is None:
            self.get_exon(gene_start)
        # for exon
        line_str = []
        for exon in self.exon:
            start, end = exon
            str_one = "\t".join([self.chrom, str(start), str(end), self.name])
            line_str.append(str_one)

        bed_str = "\n".join(line_str)
        self.exon_str=bed_str

        # re init the list for intron
        line_str = []
        for intron in self.intron:
            start, end = intron
            str_one_intron = "\t".join([self.chrom, str(start), str(end), self.name])
            line_str.append(str_one_intron)

        bed_str = "\n".join(line_str)
        self.intron_str=bed_str

    def to_bedfile(self, gene_start=None, dirpath=dirp, rewrite=False):

        exon_file=dirpath+"/"+self.name+"_exon.bed"
        intron_file=dirpath+"/"+self.name+"_intron.bed"

        has_file= os.path.isfile(exon_file) and os.path.isfile(intron_file)

        if has_file and rewrite==False:
            self.exon_file=exon_file
            self.intron_file=intron_file

        else:
            if gene_start is None:
                gene_start=0 # seems not affect the speed of bedtools
            if self.exon_str is None:
                self.to_bedstr(gene_start)

            with open(exon_file, "w") as fw:
                fw.write(self.exon_str)
            self.exon_file=exon_file

            with open(intron_file, "w") as fw:
                fw.write(self.intron_str)
            self.intron_file=intron_file

    @staticmethod
    def bedfile_cal_distance(bedfile1, bedfile2, min_length, by="ratio"):

        jaccard =wrapper_bedtools_jaccard(bedfile1, bedfile2)

        if by == "ratio":
            # intron could be 0
            if min_length == 0:
                return 0
            else:
                similar = jaccard["jaccard"]  # equals float(jaccard["intersection"])/jaccard["union-intersection"]
                return 1 - similar

        elif by == "length":
            return jaccard["union-intersection"] - jaccard["intersection"]

        elif by == "ratio_short":
            # intron could be 0
            if min_length == 0:
                return 0
            else:
                return 1 - float(jaccard["intersection"]) / min_length

        elif by == "length_short":
            return min_length - jaccard["intersection"]

    def bedfile_cal_distance_exon(self, other_bgp, gene_start=None, by="ratio"):
        """

        :param other_bgp:
        :param gene_start:
        :param by: could be "ratio", "length", "ratio_short", "length_short"
        :return: dis-similarity matrix
        """
        if self.exon_file is None:
            self.to_bedfile(gene_start)
        if other_bgp.exon_file is None:
            other_bgp.to_bedfile(gene_start)
        bed1 = self.exon_file
        bed2 = other_bgp.exon_file
        min_length = self.exonlen if self.exonlen - other_bgp.exonlen <= 0 else other_bgp.exonlen

        distance = self.bedfile_cal_distance(bed1, bed2, min_length, by)
        return distance

    def bedfile_cal_distance_intron(self, other_bgp, gene_start=None, by="ratio"):
        """
        :param other_bgp:
        :param gene_start:
        :param by: could be "ratio", "length", "ratio_short", "length_short"
        :return: dis-similarity matrix
        """
        if self.intron_file is None:
            self.to_bedfile(gene_start)
        if other_bgp.intron_file is None:
            other_bgp.to_bedfile(gene_start)

        bed1 = self.intron_file
        bed2 = other_bgp.intron_file
        min_length = self.intronlen if self.intronlen - other_bgp.intronlen <= 0 else other_bgp.intronlen

        distance = self.bedfile_cal_distance(bed1, bed2, min_length, by)
        return distance


    def bind_seq(self, seqdic):
        """
        To bind the sequence of the bigg, can used to call sl or cds frame
        :param: seqdic: A biopython seqdic that contains the sequence
        """
        self.seq= str(seqdic[self.name].seq)


    def orf_find(self, refdic):
        """
        :param: refdic: the reference genome
        """
        pass

    def write_subread(self):
        """
        subread is generated in clustering
        :return:
        """
        self.name2=",".join(self.subread)

    def get_subread_from_str(self):
        if "," in self.name2:
            self.subread=self.name2.split(",")
        elif self.name2=="none":
            self.subread=[]

    def __score_sl(self, seq1="CUCAAACUUGGGUAAUUAAACCG"):
        """
        use current best ssw aligner parameter
        write the ssw align score to score first 23 nt and write to score
        :param: default seq1 is SL1 sequence for nematodes
        """
        pass
        seq2=self.seq[:23]
        #aln = ssw_wrapper(seq1, seq2, 2, 2, 3, 1)
        #self.score = aln.score



