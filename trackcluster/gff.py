#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/8/2018 2:44 PM
# @Author  : Runsheng
# @File    : gff.py

"""
Used to parser the gff gene terms and convert to bigGenePred object
GFF is a collection of one large gff file
"""

# standard library import
from collections import namedtuple, OrderedDict
from operator import attrgetter

# set nametuple for gff class
GFF_FIELD = ['seqid', 'source', 'type', 'start', 'end',
                 'score', 'strand', 'frame', 'attributes']

GFF_record = namedtuple("gff", GFF_FIELD)


def parse_gff_line(line):
    """
    :return: gff_records list, contains 9 col nametuple
    """

    line_l = line.strip().split("\t")

    # test if the line is 9 cols
    if len(line_l) != len(GFF_FIELD):
        print(("{0} != {1}".format(len(line_l), len(GFF_FIELD))))
        pass

    # in case some gff have quotes in the cols, replace them
    else:
        data = {
            'seqid': None if line_l[0] == '.' else line_l[0].replace('"', ''),
            'source': None if line_l[1] == '.' else line_l[1].replace('"', ''),
            'type': None if line_l[2] == '.' else line_l[2].replace('"', ''),
            'start': None if line_l[3] == '.' else int(line_l[3]),
            'end': None if line_l[4] == '.' else int(line_l[4]),
            'score': None if line_l[5] == '.' else float(line_l[5]),
            'strand': None if line_l[6] == '.' else line_l[6].replace('"', ''),
            'frame': None if line_l[7] == '.' else line_l[7].replace('"', ''),
            'attributes': parse_attributes(line_l[8])
        }

        return GFF_record(**data)

def gff_record_to_str(record):
    """
    reverse function for parse_gff_line
    :param record:
    :return:
    """
    attr_str=attr_dic_to_str(record.attributes)
    line_l=list(record)
    line_l[-1]=attr_str
    line_str=["." if x is None else str(x) for x in line_l] # change the None back to .
    return "\t".join(line_str)


def parse_attributes(attr_string):
    """
    used to parse the GFF attributes field (col 9) into a dict
    example: ID=ctg7180000006767:hsp:104246:3.10.0.0;Parent=ctg7180000006767:hit:31434:3.10.0.0; \
    Target=CBG07519 47 108;Length=292;Gap=F2 M62 R1
    return: [('ID', 'ctg7180000006767:hsp:104246:3.10.0.0'), ('Parent', 'ctg7180000006767:hit:31434:3.10.0.0'), \
     ('Target', 'CBG07519 47 108'), ('Length', '292'), ('Gap', 'F2 M62 R1')]))
    """
    attr_dic = OrderedDict()

    if attr_string == ".":
        return attr_dic
    else:
        attr_dic = OrderedDict()
        for attribute in attr_string.strip().split(";"):
            if len(attribute) > 0 :
                if "=" in attribute.strip():
                    elems = attribute.strip().split('=')
                    key = elems[0]
                    value = ' '.join(elems[1:])
                    if value[0] == '"':
                        value = value[1:]
                    if value[-1] == '"':
                        value = value[0:-1]
                    attr_dic[key] = value
    return attr_dic

def attr_dic_to_str(attr_dic):
    """
    reverse function of parser_attr
    :param attr_dic:
    :return: string str for a valid gff file
    """
    if len(attr_dic)==0:
        return ""
    else:
        line_l=[]
        for k, v in attr_dic.items():
            str_one=k+"="+v
            line_l.append(str_one)
        return ";".join(line_l)



class GFF(object):
    """
    A general class to parse the gff file
    For speed issue, just read all lines into mem

    parser the gff using

    """

    def __init__(self, gfffile):

        """
        :param gfffile:
        """

        # init in file reading
        self.filename=gfffile
        self.gff_string_list = []
        self.to_list() # here use list for the flexibility instead of generator

        self.gff_records=[] # hold multiple GFF_record in a list

        # other
        self.gene_d=None
        self.transcript_to_gene=None
        self.transcript_d=None


    def to_list(self):
        """
        :return: gfflines list
        """
        with open(self.filename, "r") as f:
            lines = f.readlines()
            for line in lines:
                if line[0] != "#" and len(line) > 18: # ignore the # and small lines
                    self.gff_string_list.append(line.strip())

    def parser(self):
        if len(self.gff_string_list)==0:
            self.to_list()
        else:
            for line in self.gff_string_list:
                self.gff_records.append(parse_gff_line(line))

    def gene_format(self, indicator="ID"):
        """
        format the gff line to blocks,
        {gene name, [mRNAline, cdsline...]}

        :param indicator: usually be "ID", "Name", "sequence_name", the name used in gene line, which will be further used in exon line and/or cds line

        the gff need to be sorted and all gene/mRNA/exon/cds/UTR is together


        select the gene lines, including the
        gene, mRNA, exon,CDS, UTR
        exclude the other lines, may need to add the other lines like pseudogene,
        should try to include all keys containing an exon
        switch to bottom-up method to avoid non-exon information
        :return:
        """
        if len(self.gff_records)==0:
            self.parser()

        ####
        parent2line_dic={} # dic for {namexxx:[1,2,3...]}. collect the level 1 ID and parent information
        nontop_s=set() # set to store all lines numbers, which contains parent in attr
        id2pos={} # {namexxx:2}
        pos2id={}

        # collect all three dic
        for n, record in enumerate(self.gff_records):
            try:
                id = record.attributes[indicator]
                id2pos[id] = n
                pos2id[n]=id
            except KeyError:
                pass

            try: # record all parent
                p1=record.attributes["Parent"]
                try:
                    parent2line_dic[p1].append(n)
                except KeyError:
                    parent2line_dic[p1]=[n]
                nontop_s.add(n)
            except KeyError:
                pass

        #print(parent2line_dic)
        #print(nontop_s)
        #print(id2pos)

        # add the top key to in id2line to the gene2line
        # for ensembl gff, the top line with indicator="ID" will keep only genes
        # usually the line is one mRNA line
        gene2line={}
        for id, pos_l in parent2line_dic.items():
            try:
                if id2pos[id] not in nontop_s or len(pos_l)==0: # some gff will have missing parent line, it is wrong but it happens
                    gene2line[id]=pos_l
            except KeyError:
                pass
        #print("gene2line,",len(gene2line), "parent2line", len(parent2line_dic))
        #print_dic(gene2line)

        # reverse the line2gene
        line2gene={}
        for gene, pos_l in gene2line.items():
            for pos in pos_l:
                try:
                    record=self.gff_records[pos]
                    name=record.attributes[indicator]
                    line2gene[name]=gene
                except KeyError:
                    pass

        #add the nontop ones to top keys
        # use id2line to glue all untop line and their subline to gene2line
        # the value in gene2line need to be expaned
        for id, pos_l in parent2line_dic.items():
            try:
                if id2pos[id] in nontop_s: # not top key, has parent, need to add to gene2line
                    try:
                        gene=line2gene[id] #
                        #print("in", gene)
                        gene2line[gene].extend(pos_l)
                    except KeyError:
                        pass
            except KeyError:
                pass
        # create the gene_d
        gene_d = OrderedDict()
        for k, v in gene2line.items():
            gene_d[k]=[]
            line_gene=id2pos[k]
            gene_d[k].append(self.gff_records[line_gene])

            line_other=v
            for i in line_other:
                gene_d[k].append(self.gff_records[i])
        self.gene_d=gene_d

    def transcript_format(self, indicator="ID", keys=None):
        transcript_to_gene={} # {mrna_name: genename}
        transcript_d=OrderedDict() # {mrna_name: gff_line}

        if self.gene_d is None:
            self.gene_format(indicator=indicator)

        if keys is None: # if no chose ones, use all mrna
            keys=list(self.gene_d.keys())

        for k in keys:
            record_list=self.gene_d[k] # use record instead of gff

            # first test the type of the gene by test if CDS is in the gff_line
            type_set=set([x.type for x in record_list])
            is_protein_coding=True if "CDS" in type_set else False
            exon_junction=set()

            for record in record_list:
                if record.type == "exon":
                    transcript_name = record.attributes["Parent"]  # general
                    transcript_to_gene[transcript_name] = k

                    try:
                        transcript_d[transcript_name].append(record)
                    except KeyError:
                        transcript_d[transcript_name] = []
                        transcript_d[transcript_name].append(record)

                    exon_junction.add(record.start)
                    exon_junction.add(record.end)


        for k, v in list(transcript_d.items()):
            v.sort(key=attrgetter('start'))

        self.transcript_to_gene=transcript_to_gene
        self.transcript_d=transcript_d

    def gff_write(self, out, keys=None):
        """
        write the formatted or the original gff to a new file
        :param out:
        :param keys:
        :return:
        """

        def write_line(gff_record):
            gff_str = gff_record_to_str(gff_record)
            fw.write(gff_str)
            fw.write("\n")

        fw=open(out,"w")

        if keys is None:
           keys=list(self.gene_d.keys())

        for k in keys:
            print(k)
            record_l=self.gene_d[k]
            for record in record_l:
                write_line(record)

        fw.close()


class __WormBaseGFF(GFF):
    """
    specific parser for WormBase GFF file, include:
    un-sorted CDS and exon
    multiple parents:
    mRNA (transcript):Parent=Gene:WBxxxxx
    exon:Parent=Transcript:Y74C9A.3
    CDS: same
    intron: same
    UTR: same

    """
    pass





