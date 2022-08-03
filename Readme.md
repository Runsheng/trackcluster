# TrackCluster
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
2. python modules: pysam, pandas, numpy, biopython, tqdm
3. samtools V2.0+ , bedtools V2.24+  and minimap2 V2.24+ in your $PATH

## Recommendations
1. UCSC Kent source tree (for generating binary track)

## <a name="walkthrough"></a>Walkthrough

An walkthrough example can be found in the [ipython notebook file](https://github.com/Runsheng/trackcluster/blob/master/trackcluster_run_example.ipynb). 

The new junction mod to run trackcluster can be found in 

## Citation
Please kindly cite our paper in Genome Research if you use trackcluster in your work.

Li, R., Ren, X., Ding, Q., Bi, Y., Xie, D. and Zhao, Z., 2020. Direct full-length RNA sequencing reveals unexpected transcriptome complexity during *Caenorhabditis elegans* development. **Genome research**, 30(2), pp.287-298.