# TrackCluster
![PyPI](https://img.shields.io/pypi/v/trackcluster?color=green)

Trackcluster is an isoform calling and quantification pipeline for Nanopore direct-RNA long reads.

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Walkthrough](#walkthrough)

####
Hint: the new feather including can be found in the dev branch
- the compatibility with py3
- quick junction mode for long reads with high accuracy. 

## <a name="overview"></a>Overview
A pipeline for reference-based identification of isoform calling using Nanopore direct-RNA long reads. This pipeline was designed to use **only** long and nosisy reads to make a valid transcriptome. An indicator for the intact 5' could be very helpful to the pipeline, i.e, the splicing leader in the mRNA of nematodes. 

It is recommended to combine all samples together to generate a new transcriptome reference. After this process, the expression of isoforms in each sample can be fetched by providing an "name:sample" table. 

The output format for this pipeline is ["bigGenePred"](https://github.com/Runsheng/trackcluster/blob/master/test/bigGenePred.as). 

## <a name="requirements"></a>Requirements

1. python 3.9 (or 2.7.10+)
2. samtools V2.0+ , bedtools V2.24+  and minimap2 V2.24+ in your $PATH

## Installation
```bash
# use pip from pypi
pip install trackcluster
# or pip from source code for the latest version
git clone https://github.com/Runsheng/trackcluster.git
pip install ./trackcluster
```

## Recommendations
1. UCSC Kent source tree (for generating binary track), used only in bigg2b.py

## <a name="walkthrough"></a>Walkthrough
```bash
# test if all dependencies are installed
trackrun.py test --install

# prepare the reference annotation bed file from gff file
# tested on Ensembl, WormBase and Arapost gff
gff2bigg.py -i ensemblxxxx.gff3 -o ref.bed 
# WormBase full gff contains too many information, need to extract the lines from WormBase only
cat c_elegans.PRJNA13758.WS266.annotations.gff3 |grep WormBase > ws266.gff
gff2bigg.py -i ws266.gff -o ref.bed
# the ref.bed can be sorted to speed up the analysis
bedtools sort -i ref.bed > refs.bed # refs.bed contains the sorted, know transcripts from gff annotation

# generate the read track from minimap2 bam file
bam2bigg.py -b group1.bam -o group1.bed
bam2bigg.py -b group2.bam -o group2.bed

# merge the bed file and sort
cat group1.bed group2.bed > read.bed
bedtools sort -i read.bed > reads.bed

# Examples for running commands:
trackrun.py clusterj -s reads.bed -r refs.bed -t 40 # run in junction mode, will generate the isoform.bed
trackrun.py count -s reads.bed -r refs.bed -i isoform.bed # generate the csv file for isoform expression
# alternative for cluster
trackrun.py cluster -s reads.bed -r refs.bed -t 40 # run in exon/intron intersection mode， slower, will generate the isoform.bed

# the post analysis could include the classification of novel isoforms
trackrun.py desc --isoform isoform.bed --reference ref.bed > desc.bed  # generate the description for each novel isoform
# this part can be run directly on reads, to count the frequency of splicing events in reads, like intron_retention
trackrun.py desc --isoform reads.bed --reference ref.bed > reads_desc.bed 

```


## Citation
Please kindly cite our paper in Genome Research if you use trackcluster in your work.

Li, R., Ren, X., Ding, Q., Bi, Y., Xie, D. and Zhao, Z., 2020. Direct full-length RNA sequencing reveals unexpected transcriptome complexity during *Caenorhabditis elegans* development. **Genome research**, 30(2), pp.287-298.