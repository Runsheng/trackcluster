"""
The main flow file to process one fastq file with
"""

from trackcluster.utils import myexe, is_bin_in, get_file_prefix,del_files,parmap
from trackcluster.convert import sam_to_bigGenePred
from trackcluster.batch import process_one_junction_corrected_try
from trackcluster.tracklist import read_bigg, write_bigg

from trackcluster.pre import wrapper_bedtools_intersect2_select, tracklist_add_gene,get_gendic_bedinter,group_bigg_by_gene

#std lib
import logging
import os

#third party lib
from pysam import AlignmentFile
from tqdm import tqdm


def flow_mapping(wkdir,ref_file,fastq_file,prefix, core=16):
    """
    minimap2 mapping flow
    prefix is used in the output bam file
    :return: filename for bamfile
    """
    if is_bin_in("samtools") and is_bin_in("minimap2"):
        pass
    else:
        raise Exception("Check samtools and bedtools installion")

    os.chdir(wkdir)

    cmd_map="minimap2 -ax splice -k14 -uf -t {core} {ref} {fastq_file} | samtools view -bS -F260 -q 30 > {prefix}.bam".format(
    prefix = prefix, core = core, fastq_file = fastq_file, ref = ref_file)
    print(cmd_map)
    myexe(cmd_map)

    cmd_sam2="samtools sort -@{core} {prefix}.bam >{prefix}_s.bam".format(prefix=prefix, core=core)
    print(cmd_sam2)
    myexe(cmd_sam2)

    cmd_sam3="samtools index {prefix}_s.bam".format(prefix=prefix)
    print(cmd_sam3)
    myexe(cmd_sam3)

    del_files(["{prefix}.bam".format(prefix=prefix)])

    return "{prefix}_s.bam".format(prefix=prefix)


def flow_bamconvert(wkdir,bamfile,out,prefix,score=30):
    """
    The prepare part
    :return: write the out bed file as bigglist, return the filename
    """
    os.chdir(wkdir)
    samfile = AlignmentFile(bamfile)

    fw = open(out, "w")

    for n, record in enumerate(samfile):
        # add mapq filter to rm the secondary and supplementary mapping
        # score has been filtered in samtools view
        if record.mapq >= score:
            try:
                bigg = sam_to_bigGenePred(record, samfile)
                bigg.geneName2=prefix # add group name as geneName2
                fw.write(bigg.to_str())
                fw.write("\n")
            except ValueError:
                pass
    fw.close()
    samfile.close()
    return out


def flow_preparedir(wkdir, prefix, bigg_gff_file, bigg_nano_file):
    """
    use nanopore and gff annotation to build folders
    :param bigg_gff:
    :param bigg_nano:
    :return:  write the file for novel gene out

    write several internal files for further use:

    """
    os.chdir(wkdir)

    ### get two parts, the gene part and the novel part
    # the gene part
    outfile=prefix+"inter.bed"
    out_inter=wrapper_bedtools_intersect2_select(bigg_nano_file, bigg_gff_file, outfile=outfile,
                                                 fraction=0.2)
    read_gene=get_gendic_bedinter(outfile)
    print("read number in genes:", len(read_gene))

    bigg_nano = read_bigg(bigg_nano_file)
    bigg_new=tracklist_add_gene(bigg_nano, read_gene)

    bigg_ref= read_bigg(bigg_gff_file)

    ###
    ### create dirs for known genes
    gene_nano, bigg_novel_gene=group_bigg_by_gene(bigg_new)
    gene_anno, _ =group_bigg_by_gene(bigg_ref)

    for gene, nano_bigg in gene_nano.items():
        anno_bigg = gene_anno[gene]
        try:
            os.mkdir(gene)
        except OSError:
            pass

        anno_out = "./{gene}/{gene}_gff.bed".format(gene=gene)
        nano_out = "./{gene}/{gene}_nano.bed".format(gene=gene)

        write_bigg(anno_bigg, anno_out)
        write_bigg(nano_bigg, nano_out)
    ###
    ###

    ### the novel part, use bigg_novel_gene
    novel_file=prefix+"_novel.bed"
    write_bigg(bigg_novel_gene, novel_file)

    # write the genename to pickle
    genename_file=prefix+"_gene.txt"

    with open(genename_file, "w") as fw:
       for gene in list(gene_nano.keys()):
           fw.write(gene)
           fw.write("\n")

    return list(gene_nano.keys())


def flow_gene_clusterj(wkdir, genename_file, core=30):
    """

    :param wkdir:
    :param genename_file:
    :param core:
    :return:
    """
    os.chdir(wkdir)

    gene_l = []
    with open(genename_file, "r") as f:
        for line in f.readlines():
            gene_l.append(line.strip())

    print("###### Run gene cluster ######")
    parmap(process_one_junction_corrected_try, tqdm(gene_l), core)



def prepare_novel_dir(wkdir, novel_file):
    """
    use the novel bigg file, run merge to get the merged full bed, use each line to create a
    :return:
    """
    pass





