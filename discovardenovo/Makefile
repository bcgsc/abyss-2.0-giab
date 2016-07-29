# Number of threads
t=64

all: \
	hsapiens-scaffolds.fa \
	BESST/hsapiens-scaffolds.fa \
	abyss-scaffold/hsapiens-scaffolds.fa

# Report run time and memory usage
export SHELL=zsh -opipefail
export REPORTTIME=0
export TIMEFMT=%J  %U user %S system %P cpu %*E total %M MB

.DELETE_ON_ERROR:
.SECONDARY:

# DISCOVARdenovo

discovardenovo/$(name)/a.final/a.lines.fasta:
	time ./run-discovardenovo-assembly

$(name)-scaffolds.fa: $(name)/a.final/a.lines.fasta
	seqtk seq $< >$@

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

# ABySS

# Scaffold using abyss-scaffold

abyss-scaffold/%-6.fa: hsapiens-scaffolds.fa
	mkdir -p $(@D)
	ln -sf ../$< $@

%.dot: %.fa
	time abyss-todot $< >$@

#abyss-scaffold/%-6.sam.gz: %.bam
#	time samtools view -h -F2304 $< | samtools sort -n -OSAM | abyss-fixmate -v -h $(@D)/$*-6.hist | sort -snk3 -k4 | gzip >$@

# Align the reads to the assembled genome
abyss-scaffold/mp6k-6.sam.gz: hsapiens-scaffolds.fa \
		/projects/btl/datasets/hsapiens/giab/HG004/MPHG004-23100079/MPHG004_S3_L003_001.mp.bfc.fq.gz \
		/projects/btl/datasets/hsapiens/giab/HG004/MPHG004-23110109/MPHG004_S3_L003_001.mp.bfc.fq.gz \
		/projects/btl/datasets/hsapiens/giab/HG004/MPHG004-23100079/MPHG004_S3_L003_001.unknown.bfc.fq.gz \
		/projects/btl/datasets/hsapiens/giab/HG004/MPHG004-23110109/MPHG004_S3_L003_001.unknown.bfc.fq.gz
	time gunzip -c $(wordlist 2, 99, $^) | bwa mem -t$t -p hsapiens-scaffolds.fa - | samtools view -h -F2304 \
	| abyss-fixmate -v -h $(@D)/mp6k-6.hist | sort -snk3 -k4 | gzip >$@

%-scaffolds.fa: %-6.fa %-6.dot abyss-scaffold/mp6k-6.sam.gz
	time abyss-pe -C $(@D) -f ../intermediate.mk \
		name=hsapiens mp=mp6k \
		j=$t v=-v k=200 l=40 s=500 S=500-5000 N=15 mp6k_de=--mean mp6k_n=1 \
		scaffolds

# BESST

# Scaffold the assembled genome using BESST
BESST_output/pass1/Scaffolds_pass1.fa: hsapiens-scaffolds.fa mp6k.bam.bai
	runBESST --orientation rf -c $< -f mp6k.bam -o .

# Delete line breaks and underscores from BESST scaffolds
BESST/hsapiens-scaffolds.fa: BESST_output/pass1/Scaffolds_pass1.fa
	ln -sf BESST_output BESST
	seqtk seq BESST_output/pass1/Scaffolds_pass1.fa | tr -d _ >$@
