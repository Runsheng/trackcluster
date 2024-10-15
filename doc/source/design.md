# Trackcluster design
## History of ONT read accuracy and trackcluster
The direct-RNA sequencing from ONT used to have very low quaility, especially for the non-human samples. The raw read quality
is roughly 85% for RNA001 kit with R9.4.1 flowcell. With a 15% error rate, the read junctions are also error prone, which
makes the isoform calling and quantification using junctions very difficult. To accomendate the low quality reads, we have
designed trackcluster by comparing the intersection of exon/intron regions to determine the distance between different reads,
and try to correct the junctions after clustering all similar reads from one isoform. 

The regional intersection method (original trackcluster) worked well for RNA001/002 data in model organisms like _C. elegans_ and _Arabidopsis thaliana_, 
as the read count for each isoform is limited. The total yield for one flowcell is around 1-2Gb. And for one experiment, the 
overall yield is usually below 10Gb, and the gene expression for highly expressed genes are generally less than 50,000. 
However, the method is not fast enough as we have to calculate the intersection for every two reads who have an overlap. The
time complexity is O(n^c), which is not acceptable for genes with expression higher than 50,000 (may take 24h for computation).


With the new RNA002 kit and new basecall models, the read quality is improved to 92% for most samples. And 8% (instead of >15%) 
of error rate would allow for the junction self-correction before clustering. So we also included the junction self-correction
methods, and also the clustering methods using junctions (trackclusterj method). The time complexity for trackclusterj is roughly O(nlogn) for most 
of the cases, which is acceptable for the high expressed genes.