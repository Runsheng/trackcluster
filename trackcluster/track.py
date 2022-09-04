#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/8/2018 2:44 PM
# @Author  : Runsheng     
# @File    : track.py

# third part import
from trackcluster.utils import chr_select, reverse_complement
import os


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
    string name2; "Alternative/human readable name, |store subreads or coverage"
    string cdsStartStat; "Status of CDS start annotation (none, unknown, incomplete, or complete)"
    string cdsEndStat; "Status of CDS end annotation (none, unknown, incomplete, or complete)"
    int[blockCount] exonFrames; "Exon frame {0,1,2}, or -1 if no frame for exon"
    string type; "Transcript type"
    string geneName; "Primary identifier for gene"
    string geneName2; "Alternative/human readable gene name, |used to store sample information such as stage or group"
    string geneType; "Gene type"
    )
    """
    def __init__(self, prefix=""):
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

        # prefix was needed when used for wormbase to ucsc
        # useful to convert "I" to "chrI"
        self.prefix=prefix

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
        self.subread=set() # use to store the reads contained inside the
        self.coverage=0

        self.junction=[]

        ###
        self.seq_chro=None

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

        # ucsc compatible chr name, only use for the Worm chr as I, II....
        return self.prefix+"\t".join(str_l)

    def __str__(self):
        return self.to_str()

    def from_string(self, string_bed):
        """
        IO
        reverse function of the to_str
        """
        string_l=string_bed.strip().split("\t")
        chr_prefix=self.prefix

        self.chrom=string_l[0].replace(chr_prefix, "") # revsere the +prefix operation
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

    @staticmethod
    def get_intron(line_exon):
        line_intron=[]
        for n, pair in enumerate(line_exon):
            start, end=pair
            if n==0:
                pass
            else:
                line_intron.append((end_p, start))
            if n==len(line_exon):
                break
            start_p, end_p=pair
        return line_intron

    @staticmethod
    def get_exonlen(line_exon):
        len_exon=0
        for x,y in line_exon:
            len_exon+=(y-x)
        return len_exon

    def get_exon(self, gene_start=None):
        """
        note the 0 and 1 based sys, [0:100) is the regions
        contain the 0 , but do not contain 100
        :param gene_start:
        :return:
        """
        gene_start=0 if gene_start is None else gene_start
        offset=self.chromStart-gene_start

        line_exon=[]

        for x, y in zip(self.chromStarts, self.blockSizes):
            line_one=(x+offset, x+offset+y)
            line_exon.append(line_one)

        self.exon=line_exon
        self.exonlen=self.get_exonlen(line_exon)

        self.intron=self.get_intron(line_exon)
        self.intronlen=self.get_exonlen(self.intron)

        # debug:
        #print("exon", self.exon)

    def exon_to_block(self, gene_start=None):
        """
        reverse of get_exon
        :param gene_start:
        :return:
        """
        gene_start=0 if gene_start is None else gene_start
        offset=self.chromStart-gene_start

        start_l=[]
        block_l=[]

        for pair in self.exon:
            start, end=pair
            chrom_start=start-offset
            block_size=end-start

            start_l.append(chrom_start)
            block_l.append(block_size)

        self.chromStarts=start_l
        self.blockSizes=block_l

    def to_bedstr(self, gene_start=None):
        """
        convert the exon region to BedTool object
        :return:
        """
        if self.exon is None:
            self.get_exon(gene_start=gene_start)
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
        if self.intronlen>0:
            for intron in self.intron:
                start, end = intron
                str_one_intron = "\t".join([self.chrom, str(start), str(end), self.name])
                line_str.append(str_one_intron)
        else: #  for single exon gene,generate a line with len=1, start=end-1=chromStart
            start=self.chromStart
            end=self.chromStart+1
            str_one_intron = "\t".join([self.chrom, str(start), str(end), self.name])
            line_str.append(str_one_intron)

        bed_str = "\n".join(line_str)
        self.intron_str=bed_str

    def bind_readseq(self, seqdic):
        """
        To bind the sequence of the bigg, can used to call sl or cds frame
        :param: seqdic: A biopython seqdic that contains the sequence
        """
        self.seq= str(seqdic[self.name].seq)

    @staticmethod
    def get_cumsum(l):
        new_l = []
        cumsum = 0
        for elt in l:
            cumsum += elt
            new_l.append(cumsum)
        return new_l

    def mrna_pos_to_chro(self, mrna_pos):
        """
        Convert the position in the mrna to the position in chro
        two steps: 1, determine the exon number 2, get the pos

        ##todo: add the reverse strand code and test

        mrna poition is 0 based
        chro position is 0 based
        :param mrna_pos:
        :return: [chro(str), pos(int)]
        """
        if self.exon is None:
            self.get_exon()
        block_sum=self.get_cumsum(self.blockSizes)


        ### sanity check
        if 0<=mrna_pos<=block_sum[-1]:
            pass
        else:
            return None

        mrna_pos=mrna_pos if self.strand == "+" else (self.exonlen-mrna_pos-1)
        print(mrna_pos)

        ### for the forward mRNA
        block_sum_1=[0]+block_sum
        # define which exon is the mrna_pos in
        for i, exon_sum in enumerate(block_sum):
            if mrna_pos<=exon_sum:
                exon_num=i
                exon_offset=mrna_pos-block_sum_1[i]
                break

        chro_pos=self.chromStart+self.chromStarts[exon_num]+exon_offset

        return self.chrom, chro_pos

    def bind_chroseq(self, refdic, gap=0, intron=False):
        """
        the gap=0 and intron=False will output exon seq
        gap>0 would not work with intron=True
        gap>0 and intron=False will output exon seq seperated by "N"
        :param: refdic: the reference genome in a dict
        """
        # need self.exon, self.strand
        # need to note that the chr_select is 0 based
        # while the exon in bigg is 1 based
        if self.exon is None:
            self.get_exon()

        seq_l=[]
        for n,exon_one in enumerate(self.exon):
            chro=self.chrom
            start, end=exon_one
            _,seq=chr_select(refdic, chro, start,end)
            seq_l.append(seq.upper()) # upper case for exon
            if intron is False and n<len(self.exon)-1:
                seq_l.append(gap*"N")
            if intron and n<len(self.intron):
                intron_one=self.intron[n]
                start_i, end_i=intron_one
                _, seq_intron=chr_select(refdic, chro, start_i, end_i)
                seq_l.append(seq_intron.lower()) # lower case for intron
        seq_raw="".join(seq_l)
        if self.strand=="+":
            seq_out=seq_raw
        else:
            seq_out=reverse_complement(seq_raw)

        self.seq_chro="".join(seq_out)

    def orf_find(self, refdic):
        """
        :param: refdic: the reference genome
        """
        # need self.seq_chro
        if self.seq_chro is None:
            return None
        else:
            pass

    def find_orfs_with_trans(self, trans_table=1, min_protein_length=10):
        """
        compare the three orfs, use the longest one as new protein
        min size may need to be defined as 30aa?,so cds need to >90
        :param trans_table:
        :param min_protein_length:
        :return:
        """
        if self.seq_chro is None:
            return None
        else:
            from Bio.Seq import Seq
            seq=Seq(self.seq_chro)

            seq0=seq.translate(table=1)
            seq1=seq[1:].translate(table=1)
            seq2=seq[2:].translate(table=1)

            print((str(seq0)))
            print("======================================")
            print((str(seq1)))
            print("======================================")
            print((str(seq2)))

            answer=[seq0, seq1, seq2]


        return answer

    def is_nmd(self):
        """
        use the cutoff as: >=50nt of the last junction to define the PTC, and get the NMD
        :return:
        """
        pass

    def write_subread(self):
        """
        subread is generated in clustering, using cluster.py
        :return:
        """
        self.name2=",".join(list(self.subread))

    def write_coverage(self):
        """
        can only be run when coverage are calculated in tracklist
        :return:
        """
        self.name2=",".join(list(self.subread))+",|"+str(self.coverage)

    def get_coverage_from_str(self):
        """
        reverse the write_coverage
        :return:
        """
        if ",|" in self.name2:
            self.coverage=float(self.name2.split(",|")[-1])

    def get_subread_from_str(self):
        if "," in self.name2:
            self.subread=set(self.name2.split(",|")[0].split(","))
        elif self.name2=="none":
            self.subread=set()

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

    def get_junction(self):
        if self.exon is None:
            self.get_exon()

        junction_l=[]
        for exon in self.exon:
            start,end=exon
            junction_l.append(start)
            junction_l.append(end)

        # Note: the junction is sorted from the first to the last
        # Note as corrd
        #junction_lf= junction_l[1:-1] if self.strand=="+" else (junction_l[1:-1])[::-1]
        junction_lf= junction_l if self.strand=="+" else list(reversed(junction_l))

        self.junction=junction_lf[1:-1] # rm the chromStart and chromEnd

    def write_junction_to_exon(self):
        # used only for clustercj
        # sanity check
        if len(self.junction)%2!=0:
            raise SyntaxError("Junction number is not valid, need to be 2n")
            return None

        line_exon=[]
        start=self.chromStart
        end=self.chromEnd

        # remember to sort junction back to corrd
        junction=sorted(self.junction)

        for n,i in enumerate(junction):
            if n==0:
                pair=(start, i)
                line_exon.append(pair)

            elif 0<n<len(junction)-1:
                if n % 2==1:
                    pair=(i, junction[n+1])
                    line_exon.append(pair)

            elif n==len(junction)-1: # last
                pair=(i, end)
                line_exon.append(pair)

        self.exon=line_exon
        self.exonlen=self.get_exonlen(line_exon)
        self.intron=self.get_intron(line_exon)
        self.intronlen=self.get_exonlen(self.intron)














