# ABySS 2.0 assemblies of HG004

This directory contains the assemblies of HG004 from *ABySS 2.0: Resource-Efficient Assembly of Large Genomes using a Bloom Filter*. For the preprint, see <http://dx.doi.org/10.1101/068338>.

We downloaded the Illumina WGS 2x250 bp paired-end sequencing data and the Illumina 6 kbp mate-pair sequencing data of the Ashkenazi mother (HG004). We removed adapters from the mate-pair reads using NxTrim 0.4.0. We corrected sequencing errors in the reads using the tool BFC. The data is downloaded from <https://github.com/genome-in-a-bottle/giab_data_indexes>.

We assembled this data using

+ ABySS 1.9.0
+ ABySS 2.0.0 (branch bloom-abyss-preview)
+ BCALM2 2.0.0
+ DISCOVARdenovo 52488
+ MEGAHIT 1.0.6
+ Minia 3.0.0-alpha1
+ SGA 0.10.14
+ SOAPdenovo 2.04

For the assembly scripts, see the `Makefile` in <https://github.com/sjackman/giab>.
