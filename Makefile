# Assemble Genome in a Bottle
# Written by Shaun Jackman <sjackman@gmail.com>

# Name of the species and assemblies
name=hsapiens

# Number of threads
t=32

# Reference genome
ref=GRCh38
ref_fa=/projects/btl/reference_genomes/H_sapiens/GRCh38/GCA_000001405.15_GRCh38_genomic.chr-only.fa
ref_gff=/projects/btl/reference_genomes/H_sapiens/GRCh38/Homo_sapiens.GRCh38.84.chr.gff

# Parallelize gzip
gzip = pigz -p$t

# Report run time and memory usage
export SHELL=zsh -opipefail
export REPORTTIME=1
export TIMEFMT=%J  %U user %S system %P cpu %*E total %M MB

.DELETE_ON_ERROR:
.SECONDARY:

all: curl kmerstream nxtrim bfc bfc_kmerstream abyss_ref \
	bcalm discovardenovo sga soapdenovo \
	quast rmarkdown

curl: \
	AshkenazimTrio/sequence.index.AJtrio_Illumina300X_wgs_07292015.tsv \
	AshkenazimTrio/sequence.index.AJtrio_Illumina_2x250bps_02192016.tsv \
	AshkenazimTrio/sequence.index.AJtrio_Illumina_6kb_matepair_wgs_08032015.tsv \
	AshkenazimTrio/sequence.index.AJtrio_NIST_Stanford_Moleculo_125bps_08042015.tsv \
	AshkenazimTrio/alignment.index.AJtrio_10XGenomics_bwamem_GRCh37_08142015.tsv \
	HG004/sequence.index.AJtrio_Illumina_2x250bps_02192016.curl \
	HG004/MPHG004-23100079/sequence.index.AJtrio_Illumina_6kb_matepair_wgs_08032015.curl \
	HG004/MPHG004-23110109/sequence.index.AJtrio_Illumina_6kb_matepair_wgs_08032015.curl

fastqc: \
	HG004/sequence.index.AJtrio_Illumina_2x250bps_02192016.fastqc/log \
	HG004/MPHG004-23100079/sequence.index.AJtrio_Illumina_6kb_matepair_wgs_08032015.fastqc/log \
	HG004/MPHG004-23110109/sequence.index.AJtrio_Illumina_6kb_matepair_wgs_08032015.fastqc/log

kmerstream: \
	HG004/mp6k.mp.kmerstream.tsv \
	HG004/sequence.index.AJtrio_Illumina_2x250bps_02192016.kmerstream.tsv

nxtrim: \
	HG004/MPHG004-23100079/MPHG004_S3_L003_001.mp.fastq.gz \
	HG004/MPHG004-23110109/MPHG004_S3_L003_001.mp.fastq.gz

bfc: \
	HG004/sequence.index.AJtrio_Illumina_2x250bps_02192016.bfc \
	HG004/mp6k.mp.bfc.log \
	HG004/mp6k.unknown.bfc.log \
	HG004/sequence.index.AJtrio_Illumina_2x250bps_02192016.bfc.log

bfc_kmerstream: \
	HG004/mp6k.mp.bfc.kmerstream.tsv \
	HG004/sequence.index.AJtrio_Illumina_2x250bps_02192016.bfc.kmerstream.tsv

.PHONY: bcalm discovardenovo sga soapdenovo

bcalm: \
	bcalm/hsapiens-unitigs.stats.tsv \
	bcalm/GRCh38_hsapiens-unitigs.samtobreak.tsv \

discovardenovo: \
	discovardenovo/hsapiens-scaffolds.stats.tsv \
	discovardenovo/GRCh38_hsapiens-scaftigs.samtobreak.tsv \

sga: \
	sga/hsapiens-contigs.stats.tsv \
	sga/GRCh38_hsapiens-contigs.samtobreak.tsv \

soapdenovo: \
	soapdenovo/hsapiens-scaffolds.stats.tsv \
	soapdenovo/GRCh38_hsapiens-scaftigs.samtobreak.tsv \

abyss_ref: \
	$(ref)-k128/$(name)-1.fa

bionano: \
	discovardenovo/BESST/bionano/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
	discovardenovo/links/bionano/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
	discovardenovo/bionano/GRCh38_hsapiens-scaftigs.samtobreak.tsv

quast: \
	abyss/k96/$(name)-1.quast/transposed_report.tsv \
	abyss/k128/$(name)-1.quast/transposed_report.tsv \
	abyss/k160/$(name)-1.quast/transposed_report.tsv \
	abyss/k192/$(name)-1.quast/transposed_report.tsv

rmarkdown: \
	assembly-stats.html

# Install dependencies

brew_deps=abyss allpaths-lg bcalm bfc bwa discovardenovo fastqc jq kmerstream links-scaffolder masurca miller minia nxtrim pigz samtools seqtk sga soapdenovo
pip_deps=besst

deps:
	@echo $(brew_deps) $(pip_deps)

install-deps:
	brew install $(brew_deps)
	pip install $(pip_deps)

versions.json:
	brew info --json=v1 $(brew_deps) >$@

versions.tsv: %.tsv: %.json
	jq -r '["Name", "Version"],(.[]|[.name,.versions.stable])|@tsv' $< >$@

# Download the data

AshkenazimTrio/%.tsv:
	mkdir -p $(@D)
	curl https://raw.githubusercontent.com/genome-in-a-bottle/giab_data_indexes/master/AshkenazimTrio/$* >$@

HG004/%.url: AshkenazimTrio/%.tsv
	mkdir -p $(@D)
	awk '$$5 == "HG004" { print $$1 "\n" $$3 }' $< >$@

HG004/MPHG004-23100079/%.url: AshkenazimTrio/%.tsv
	mkdir -p $(@D)
	awk '/MPHG004-23100079/ { print $$1 "\n" $$3 }' $< >$@

HG004/MPHG004-23110109/%.url: AshkenazimTrio/%.tsv
	mkdir -p $(@D)
	awk '/MPHG004-23110109/ { print $$1 "\n" $$3 }' $< >$@

HG004/%.md5: AshkenazimTrio/%.tsv
	awk '$$5 == "HG004" { sub(".*/", "", $$1); sub(".*/", "", $$3); print $$2 "\t" $$1 "\n" $$4 "\t" $$3 }' $< >$@

HG004/MPHG004-23100079/%.md5: AshkenazimTrio/%.tsv
	awk '/MPHG004-23100079/ { sub(".*/", "", $$1); sub(".*/", "", $$3); print $$2 "\t" $$1 "\n" $$4 "\t" $$3 }' $< >$@

HG004/MPHG004-23110109/%.md5: AshkenazimTrio/%.tsv
	awk '/MPHG004-23110109/ { sub(".*/", "", $$1); sub(".*/", "", $$3); print $$2 "\t" $$1 "\n" $$4 "\t" $$3 }' $< >$@

%.curl: %.url %.md5
	(cd $(@D) && xargs -n1 curl -O'#') <$< 2>&1 | tee $@
	cd $(@D) && md5sum -c $(*F).md5
	cd $(@D) && for i in `cut -f2 $(*F).md5`; do ln -s $$i `echo $$i | sed -Ee 's/(.*)_R([12])(_.*).fastq.gz/\1\3_\2.fq.gz/'`; done

%.path: %.md5
	cut -f2 $< >$@

%.realpath: %.path
	(cd $(<D) && realpath `<$(<F)`) >$@

# FastQC

# Check read quality
%.fastqc/log: %.realpath
	mkdir -p $(@D)
	fastqc -t $t -o $(@D) `<$<` 2>&1 | tee $@

# KmerStream

# Count k-mers
%.kmerstream.orig: %.realpath
	KmerStream -t$t -s1 --tsv -k16,32,48,64,80,96,112,128,144,160,176,192,208,224,240 -o$@ `<$<`

# Estimate additional parameters.
%.kmerstream.tsv: %.kmerstream.orig
	awk '$$3 > 0' $< | KmerStreamEstimate.py /dev/stdin >$@

# NxTrim

# Trim mate-paire reads
%.mp.fastq.gz: %_1.fq.gz %_2.fq.gz
	nxtrim --norc --joinreads --preserve-mp -1 $*_1.fq.gz -2 $*_2.fq.gz -O $*

# BFC

bfc=HG004/sequence.index.AJtrio_Illumina_2x250bps_02192016.bfc

# Build the Bloom filter
$(bfc): %.bfc: %.realpath
	gunzip -c `<$<` | bfc -t$t -s3G -d $*.bfc /dev/stdin /dev/null

# Correct a single FASTQ file
%.bfc.fq.gz: %.fastq.gz $(bfc)
	bfc -t$t -s3G -r $(bfc) /dev/null $< | tr '\t' ' ' | $(gzip) >$@

# Correct the reads
%.bfc.log: %.path $(bfc)
	$(MAKE) `sed 's/^/$(<D)\//;s/.fastq.gz/.bfc.fq.gz/' $<` 2>&1 | tee $@

%.bfc.path: %.path
	sed 's/.fastq.gz/.bfc.fq.gz/' $< >$@

# seqtk

%.interleave.path: %.path
	cd $(<D) && for i in `<$(<F)`; do ln -s $$i `echo $$i | sed -Ee 's/(.*)_R([12])(_.*).bfc.fq.gz/\1\3_\2.bfc.fq.gz/'`; done
	sed -n 's/_R1//p' $< >$@

# Interleave paired-end reads and convert lower case nucleotides to upper case
%.bfc.fq.gz: %_1.bfc.fq.gz %_2.bfc.fq.gz
	seqtk mergepe $^ | seqtk seq -U | $(gzip) >$@

%.fa.path: %.path
	sed 's/fq.gz$$/fa.gz/' $< >$@

# Convert FASTQ to FASTA and convert lower case nucleotides to upper case
%.fa.gz: %.fq.gz
	seqtk seq -AU $< | $(gzip) >$@

# abyss-mergepairs

%.merged.path: %.path
	sed -n 's/bfc.fq.gz/bfc.merged.fq.gz/;s/_R1//p' $< >$@

# Merge overlapping paired-end reads using abyss-mergepairs
%.bfc.merged.fq.gz: %_1.bfc.fq.gz %_2.bfc.fq.gz
	abyss-mergepairs -v -p0.9 -m10 -q10 -o $*.bfc $^ >$*_merged.tsv
	mv $*.bfc_merged.fastq $*.bfc.merged.fq
	mv $*.bfc_reads_1.fastq $*.bfc.unmerged_1.fq
	mv $*.bfc_reads_2.fastq $*.bfc.unmerged_2.fq
	$(gzip) $*.bfc.merged.fq
	$(gzip) $*.bfc.unmerged_1.fq
	$(gzip) $*.bfc.unmerged_2.fq

# Miscellaneous

HG004/mp6k.mp.path:
	printf "MPHG004-23100079/MPHG004_S3_L003_001.mp.fastq.gz\nMPHG004-23110109/MPHG004_S3_L003_001.mp.fastq.gz\n" >$@

HG004/mp6k.unknown.path:
	printf "MPHG004-23100079/MPHG004_S3_L003_001.unknown.fastq.gz\nMPHG004-23110109/MPHG004_S3_L003_001.unknown.fastq.gz\n" >$@

# ABySS

# Assemble the reference genome
$(ref)-k%/$(name)-1.fa: $(ref_fa)
	mkdir -p $(ref)-k$*
	ABYSS -v -k$* -e0 -t0 -c0 $< -o $@ -s $(ref)-k$*/$(name)-bubbles.fa 2>&1 |tee $@.log

# Convert scaffolds to scaftigs
%-scaftigs.fa: %-scaffolds.fa
	seqtk seq $< | tr _ '~' | abyss-fatoagp -f $@ >$@.agp

# Calculate assembly contiguity metrics with abyss-fac
%.stats.tsv: %.fa
	abyss-fac -e3088269832 -t500 $< >$@

# Calculate assembly contiguity and correctness metrics
%.samtobreak.txt: %.sam
	(echo '==> $< <=='; bin/abyss-samtobreak-G3088269832 -l500 $<) >$@

# Convert samtobreak.txt to TSV
%.samtobreak.tsv: %.samtobreak.txt
	abyss-samtobreak-to-tsv $< >$@

# BCALM

bcalm/hsapiens-unitigs.fa:
	make -C bcalm

# DISCOVARdenovo

discovardenovo/hsapiens-scaffolds.fa:
	make -C discovardenovo

# SGA

sga/hsapiens-contigs.fa:
	make -C sga

# SOAPdenovo

soapdenovo/hsapiens-scaffolds.fa:
	make -C soapdenovo

# BioNano Genomics

bionano_prefix=EXP_REFINEFINAL1_q_bppAdjust_cmap_hsapiens-scaffolds_fa_NGScontigs_HYBRID_SCAFFOLD

# HybridScaffold
%/bionano/hybrid_scaffolds/$(bionano_prefix)_NOT_SCAFFOLDED.fasta %/bionano/hybrid_scaffolds/$(bionano_prefix).fasta: %/hsapiens-scaffolds.fa
	/usr/bin/perl /gsc/btl/linuxbrew/Cellar/iryssolve/2.1.5063/scripts/HybridScaffold/hybridScaffold.pl \
		-n $< \
		-b bionano/aggressive-B2-N2/align0/EXP_REFINEFINAL1_q.cmap \
		-c bionano/hybridScaffold_config_aggressive.xml \
		-B2 -N2 \
		-o $*/bionano \
		-r /gsc/btl/linuxbrew/Cellar/iryssolve/2.1.5063/bin/RefAligner

# Delete undescores from sequence identifiers for abyss-samtobreak
%/hsapiens-scaffolds.fa: %/hybrid_scaffolds/$(bionano_prefix)_NOT_SCAFFOLDED.fasta %/hybrid_scaffolds/$(bionano_prefix).fasta
	cat $^ | tr -d _ >$@

# BWA

# Align scaffolds to the reference using BWA-MEM
$(ref)_%.sam: %.fa
	bwa mem -xintractg -t$t $(ref_fa) $< >$@

# Samtools

# Sort a SAM file and output a BAM file
%.sort.bam: %.sam
	samtools sort -@$t -o $@ $<

# Index a BAM file
%.sort.bam.bai: %.sort.bam
	samtools index $<

# Assembly stats

assembly-stats.tsv: \
		abyss/k144/hsapiens-scaftigs.stats.tsv \
		abyss/k144/hsapiens-scaffolds.stats.tsv \
		bcalm/hsapiens-unitigs.stats.tsv \
		discovardenovo/hsapiens-scaftigs.stats.tsv \
		discovardenovo/hsapiens-scaffolds.stats.tsv \
		sga/hsapiens-contigs.stats.tsv \
		soapdenovo/hsapiens-scaftigs.stats.tsv \
		soapdenovo/hsapiens-scaffolds.stats.tsv
	mlr --tsvlite cat $^ >$@

samtobreak.tsv: \
		abyss/k144/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		bcalm/GRCh38_hsapiens-unitigs.samtobreak.tsv \
		discovardenovo/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		sga/GRCh38_hsapiens-contigs.samtobreak.tsv \
		soapdenovo/GRCh38_hsapiens-scaftigs.samtobreak.tsv
	mlr --tsvlite cat $^ >$@

%.tsv.md: %.tsv
	mlr --itsvlite --omd cat $< >$@

# QUAST

# Analayze the assemblies using QUAST
%.quast/transposed_report.tsv: %.fa
	quast.py -t$t -eL -o $*.quast -R $(ref_fa) -G $(ref_gff) $<

# R Markdown

# Render HTML from RMarkdown
%.html: %.rmd %.tsv
	Rscript -e 'rmarkdown::render("$<", "html_document")'
