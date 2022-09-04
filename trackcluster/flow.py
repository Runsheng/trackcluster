"""
The main flow file to process one fastq file with
"""

from trackcluster.utils import myexe, is_bin_in, get_file_prefix,del_files,parmap, list2file, file2list
from trackcluster.convert import sam_to_bigGenePred
from trackcluster.batch import process_one_junction_corrected_try, process_one_subsample_try
from trackcluster.tracklist import read_bigg, write_bigg, cat_bed, bigg_count_write_native, is_a_read

from trackcluster.pre import wrapper_bedtools_intersect2_select, tracklist_add_gene,get_gendic_bedinter,group_bigg_by_gene, wrapper_bedtools_merge, mergedbed2bigg, wrapper_bedtools_subtract
from trackcluster.post import flow_desc, flow_class4



#std lib
import logging
import os
import functools

#third party lib
from pysam import AlignmentFile
from tqdm import tqdm

############# shared flows for both cluster and clusterj
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

    cmd_sam2="samtools sort -/@{core} {prefix}.bam >{prefix}_s.bam".format(prefix=prefix, core=core)
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


def flow_preparedir(wkdir, prefix, bigg_gff_file, bigg_nano_file, genename_file="gene.txt",
                    f1=0.01, f2=0.05):
    """
    use nanopore and gff annotation to build folders
    :param bigg_gff:
    :param bigg_nano:
    :return:  write the file for novel gene out, return the name list for internal check

    write several internal files for further use:

    """
    os.chdir(wkdir)

    ### get two parts, the gene part and the novel part
    # the gene part
    outfile=prefix+"inter.bed"
    # write the outfile to disk
    wrapper_bedtools_intersect2_select(bigg_nano_file, bigg_gff_file, outfile=outfile,
                                                 fraction_bed1=f1, fraction_bed2=f2)
    read_gene=get_gendic_bedinter(outfile)
    print("read number in genes:", len(read_gene))

    bigg_nano = read_bigg(bigg_nano_file)
    bigg_new=tracklist_add_gene(bigg_nano, read_gene)

    bigg_ref= read_bigg(bigg_gff_file)

    ###
    ### create dirs for known genes
    gene_nano, bigg_novel_gene=group_bigg_by_gene(bigg_new)
    print("reads do not belong to gene regions: ", len(bigg_novel_gene))
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

    # write the genename to file
    genename_l=list(gene_nano.keys())

    list2file(genename_l, genename_file)

    # cleanup
    del_files([outfile])

    return genename_l


def flow_prepare_novel_dir(wkdir, prefix, novel_file, bigg_anno, genename_file="novelname.txt"):
    """
    use the novel bigg file, run merge to get the merged full bed, use each line to create a
    run the
    :return:
    """
    os.chdir(wkdir)

    # wrapper use file to transfer data
    mergedbed=prefix+"_merge.bed"
    out=wrapper_bedtools_merge(novel_file,mergedbed) # out=mergedbed
    print(out)

    bigg_regionmark=prefix+"_regionmark.bed"
    bigg_l=mergedbed2bigg(mergedbed, count_cutoff=5)
    write_bigg(bigg_l, bigg_regionmark) # write bigg_regionmark bed file

    bigg_regionmark_f = prefix + "_regionmarkf.bed"
    wrapper_bedtools_subtract(bigg_regionmark, bigg_anno, bigg_regionmark_f) # write bigg_regionmark_f file

    # create all dir and write the index file to novelname.txt
    flow_preparedir(wkdir, prefix, bigg_regionmark, novel_file, genename_file=genename_file) # write novelname.txt

    return genename_file


def flow_count(wkdir, prefix, nano_bed, isoform_bed, gff_bed):
    """
    make the count for group/sample and isoform/gene
    write the output to gene.csv and isoform.csv
    :param bigg_isoform: the bigg with subreads in the subread col
    :return: write the gene/isoform level expression
    """
    os.chdir(wkdir)

    # gff file is used to filter the reference from the read count
    bigg_gff=read_bigg(gff_bed)
    refname_set=set([bigg.name for bigg in bigg_gff])

    # the nano file is read again to get the group for each read
    bigg_nano=read_bigg(nano_bed)
    name2group={bigg.name:bigg.geneName2 for bigg in bigg_nano}
    # levels of groups, used for expression list
    groups=list(set(list(name2group.values())))
    groupvalue2pos={x:n for n,x in enumerate(groups)}
    print("groups are: ", groups)

    # isoform counting for each isoform (can use the >5 exp isoforms)
    bigg_isoform=read_bigg(isoform_bed)

    # use a list for store [(geneName, name, coverage, group1, group2...),()...]
    expression_list = []  # store

    for bigg in tqdm(bigg_isoform):
        bigg.get_coverage_from_str()
        bigg.get_subread_from_str()

        coverage = bigg.coverage
        line_prefix=[bigg.geneName, bigg.name, coverage]
        line_count=[0]*len(groups)


        for read in bigg.subread:
            if is_a_read(read,refname_set):
                try:
                    group = name2group[read]
                    pos=groupvalue2pos[group]
                    line_count[pos]+=1
                except KeyError:
                    pass # the ref isoform should be skipped

        # print g1,g2,g3
        gsum = float(sum(line_count))
        if gsum < 0.0001:
            line_one=line_prefix+[0]*len(groups)
        else:
            # adjust the count of each isoform, considering the dups in subreads
            line_count_adj=[x/gsum*coverage for x in line_count]
            line_one=line_prefix+line_count_adj
        expression_list.append(line_one)

    out=prefix+"_exp.csv"
    line_head = ["gene", "isoform", "coverageall"] + groups
    with open(out, "w") as fw:
        fw.write(",".join(line_head))
        fw.write("\n")
        for line in expression_list:
            line = [str(x) for x in line]
            fw.write(",".join(list(line)))
            fw.write("\n")

########################################################################################################


########################################################################################################
###################### flow functions for clusterj
def flow_key_clusterj(wkdir, genename_file, core=30, batchsize=2000):
    """
    run clusterj in all prepared folders, folder name from the genename_file
    can be used in both gene and landmark novel gene runs
    :param wkdir:
    :param genename_file:the namelist key file for the folders
    :param core:
    :return:
    """
    os.chdir(wkdir)

    gene_l=file2list(genename_file)
    process_one=functools.partial(process_one_junction_corrected_try, batchsize=batchsize)

    print("###### Run junction cluster ######")
    parmap(process_one, tqdm(gene_l), core)


def flow_clusterj_all_gene_novel(wkdir, prefix,nano_bed, gff_bed, core=30,
                                 f1=0.01, f2=0.01, count_cutoff=5, batchsize=2000):
    """
    wkdir and prefix will be passed from external bash wrappers
    :param wkdir:
    :param prefix:
    :param nano_bed:
    :param gff_bed:
    :param core:
    :param f1:
    :param f2:
    :return: run all cluster job
    write all inter and out file to disk
    """

    # internal file name
    novel_bed = prefix + "_novel.bed"
    mergefile = prefix + "_merge.bed"
    genename_file = prefix + "_gene.txt"
    bigg_regionmark = prefix + "_regionmark.bed"
    bigg_regionmark_f = prefix + "_regionmarkf.bed"
    novelname_file = prefix + "_novelname.txt"

    # output file
    bigg_isoform_file=prefix+"_isoform.bed"
    bigg_isoform_cov5_file=prefix+"_cov{}_isoform.bed".format(count_cutoff)

    os.chdir(wkdir)

    # step1
    print("Step1, Running prepare dir")
    # the novel bigg file is prefix_novel.bed
    # the gene name file for further use is prefix_gene.txt
    gene_l = flow_preparedir(wkdir, prefix, bigg_gff_file=gff_bed, bigg_nano_file=nano_bed,
                             genename_file=genename_file,
                             f1=f1, f2=f2)
    print("Gene name format example: ", gene_l[0])

    # step2 use genename to run gene cluster
    print("Step2, Running cluster junction")
    flow_key_clusterj(wkdir, genename_file, core=core, batchsize=batchsize)

    # combine the bed file and write the new isoforms out
    bigg_isoform=cat_bed("**/*_simple_coveragej.bed") # use ** for all file in the wkdir
    write_bigg(bigg_isoform,bigg_isoform_file)

    bigg_isoform_cov5=[]
    for bigg in bigg_isoform:
        bigg.get_coverage_from_str()
        if bigg.coverage>=count_cutoff:
            bigg_isoform_cov5.append(bigg)
    write_bigg(bigg_isoform_cov5,bigg_isoform_cov5_file)

    # step 3, prepare and run the novel bed
    print("Step3, Running novel gene finding and clustering")
    wrapper_bedtools_merge(novel_bed, mergefile)  # write mergefile
    bigg_l = mergedbed2bigg(mergefile, count_cutoff=count_cutoff)
    write_bigg(bigg_l, bigg_regionmark)  # write bigg_regionmark

    # 421 for the new genes,
    # novel should be subtract from the newly annotated isoforms to avoid the 3' short novel isoforms
    # only two left, need to adjust again
    # try to use the count >5 isoforms to avoid long intron isoform impact
    wrapper_bedtools_subtract(bigg_regionmark, bigg_isoform_cov5_file, bigg_regionmark_f, f1=f1, f2=f2)  # write bigg_regionmark_f

    # the substract will change the original novel file
    flow_preparedir(wkdir, prefix, bigg_regionmark_f, novel_bed, genename_file=novelname_file)  # write genename_file
    flow_key_clusterj(wkdir, novelname_file, core=core, batchsize=batchsize)

    return 1
########################################################################################################



########################################################################################################
###################### flow functions for cluster， which is similar but has difference in wrapper

def flow_key_cluster(wkdir, genename_file, core=30, batchsize=2000, intronweight=0.5,
                     cutoff1=0.05, cutoff2=0.005, scorecutoff=11):
    """
    run clusterj in all prepared folders, folder name from the genename_file
    can be used in both gene and landmark novel gene runs
    :param wkdir:
    :param genename_file:the namelist key file for the folders
    :param core:
    :param cutoff1 and 2: cutoff used for round1 cluster and round2 cluster
    :return:
    """
    os.chdir(wkdir)
    gene_l=file2list(genename_file)
    process_one=functools.partial(process_one_subsample_try,
                                  batchsize=batchsize, intronweight=intronweight,
                                  cutoff1=cutoff1, cutoff2=cutoff2, scorecutoff=scorecutoff)

    print("###### Run intersection cluster ######")
    parmap(process_one, tqdm(gene_l), core)



def flow_cluster_all_gene_novel(wkdir, prefix,nano_bed, gff_bed, core=30,f1=0.01, f2=0.05, count_cutoff=5,
                                batchsize=2000, intronweight=0.5,cutoff1=0.05, cutoff2=0.005, scorecutoff=11):
    """
    wkdir and prefix will be passed from external bash wrappers
    the original intersection function will be used
    :return: run all cluster job
    write all inter and out file to disk
    """

    # internal file name
    novel_bed = prefix + "_novel.bed"
    mergefile = prefix + "_merge.bed"
    genename_file = prefix + "_gene.txt"
    bigg_regionmark = prefix + "_regionmark.bed"
    bigg_regionmark_f = prefix + "_regionmarkf.bed"
    novelname_file = prefix + "_novelname.txt"

    # output file
    bigg_isoform_file=prefix+"_isoform.bed"
    bigg_isoform_cov5_file=prefix+"_cov{}_isoform.bed".format(count_cutoff)

    os.chdir(wkdir)

    # step1
    print("Step1, Running prepare dir")
    # the novel bigg file is prefix_novel.bed
    # the gene name file for further use is prefix_gene.txt
    gene_l = flow_preparedir(wkdir, prefix, bigg_gff_file=gff_bed, bigg_nano_file=nano_bed,
                             genename_file=genename_file,
                             f1=f1, f2=f2)
    print("Gene name format example: ", gene_l[0])

    # step2 use genename to run gene cluster, the function is different
    # modify: replace the function with new parameters

    flow_key_cluster_use=functools.partial(flow_key_cluster, core=core, batchsize=batchsize,
                                           intronweight=intronweight,
                                           cutoff1=cutoff1, cutoff2=cutoff2, scorecutoff=scorecutoff)

    print("Step2, Running cluster junction")
    flow_key_cluster_use(wkdir, genename_file)

    # combine the bed file and write the new isoforms out
    # modify: the default output surfix is different, so change to simple_coverage.bed
    bigg_isoform=cat_bed("**/*_simple_coverage.bed") # use ** for all file in the wkdir
    write_bigg(bigg_isoform,bigg_isoform_file)

    bigg_isoform_cov5=[]
    for bigg in bigg_isoform:
        bigg.get_coverage_from_str()
        if bigg.coverage>=count_cutoff:
            bigg_isoform_cov5.append(bigg)
    write_bigg(bigg_isoform_cov5,bigg_isoform_cov5_file)

    # step 3, prepare and run the novel bed
    print("Step3, Running novel gene finding and clustering")
    wrapper_bedtools_merge(novel_bed, mergefile)  # write mergefile
    bigg_l = mergedbed2bigg(mergefile, count_cutoff=count_cutoff)
    write_bigg(bigg_l, bigg_regionmark)  # write bigg_regionmark

    # 421 for the new genes,
    # novel should be subtract from the newly annotated isoforms to avoid the 3' short novel isoforms
    # only two left, need to adjust again
    # try to use the count >5 isoforms to avoid long intron isoform impact
    wrapper_bedtools_subtract(bigg_regionmark, bigg_isoform_cov5_file, bigg_regionmark_f, f1=f1, f2=f2)  # write bigg_regionmark_f

    # the substract will change the original novel file
    flow_preparedir(wkdir, prefix, bigg_regionmark_f, novel_bed, genename_file=novelname_file)  # write genename_file
    flow_key_cluster_use(wkdir, novelname_file)

    return 1

########################################################################################################
########################################################################################################

def flow_annotation(wkdir, ):
    pass


############
def flow_map_convert_clusterj_count(wkdir, prefix, ref_fasta, fastq_l, nano_bed, gff_bed, core=30,
                                 f1=0.01, f2=0.01, count_cutoff=5, batchsize=2000):
    """
    The overall full run pipeline for impatient people
    :param wkdir:
    :param prefix:
    :param ref_fasta:
    :param fastq_l:
    :param nano_bed:
    :param gff_bed:
    :param core:
    :param f1:
    :param f2:
    :param count_cutoff:
    :param batchsize:
    :return:
    """
    pass
