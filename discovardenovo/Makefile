# Number of threads
t=12

all: BESST/hsapiens-scaffolds.fa

.DELETE_ON_ERROR:
.SECONDARY:

# BWA

# Index the FASTA file file
%.fa.bwt: %.fa
	bwa index $<

# Align the reads to the assembled genome
mp6k.bam: hsapiens-scaffolds.fa.bwt \
		/projects/btl/datasets/hsapiens/giab/HG004/MPHG004-23100079/MPHG004_S3_L003_001.mp.bfc.fq.gz \
		/projects/btl/datasets/hsapiens/giab/HG004/MPHG004-23110109/MPHG004_S3_L003_001.mp.bfc.fq.gz \
		/projects/btl/datasets/hsapiens/giab/HG004/MPHG004-23100079/MPHG004_S3_L003_001.unknown.bfc.fq.gz \
		/projects/btl/datasets/hsapiens/giab/HG004/MPHG004-23110109/MPHG004_S3_L003_001.unknown.bfc.fq.gz
	gunzip -c $(wordlist 2, 99, $^) | bwa mem -t$t -p hsapiens-scaffolds.fa - | samtools sort -@$t -o $@

# Index the BAM file
%.bam.bai: %.bam
	samtools index $<

# BESST

# Scaffold the assembled genome using BESST
BESST_output/pass1/Scaffolds_pass1.fa: hsapiens-scaffolds.fa mp6k.bam.bai
	runBESST --orientation rf -c $< -f mp6k.bam -o .

# Delete line breaks and underscores from BESST scaffolds
BESST/hsapiens-scaffolds.fa: BESST_output/pass1/Scaffolds_pass1.fa
	ln -sf BESST_output BESST
	seqtk seq BESST_output/pass1/Scaffolds_pass1.fa | tr -d _ >$@
