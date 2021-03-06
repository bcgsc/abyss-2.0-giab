# Assemble Genome in a Bottle
# Written by Shaun Jackman <sjackman@gmail.com>

# Name of the species and assemblies
name=hsapiens

# Number of threads
t=64

# Reference genome
ref=GRCh38
ref_fa=/projects/btl/reference_genomes/H_sapiens/GRCh38/GCA_000001405.15_GRCh38_genomic.chr-only.fa
ref_gff=/projects/btl/reference_genomes/H_sapiens/GRCh38/Homo_sapiens.GRCh38.86.chr.gff3

# Parallelize gzip
gzip = pigz -p$t

# Path to ABySS executables
abyss_bin=/home/benv/arch/genesis/abyss-1.9.0/k256/bin

# Size of the reference genome with Ns
GwithN=3088269832

# Size of the reference genome without Ns
GwithoutN=2937639113

# Report run time and memory usage
export SHELL=zsh -opipefail
export REPORTTIME=1
export TIMEFMT=%J  %U user %S system %P cpu %*E total %M MB

.DELETE_ON_ERROR:
.SECONDARY:

all: curl kmerstream nxtrim bfc bfc_kmerstream abyss_ref \
	abyss bcalm discovardenovo sga soapdenovo \
	rmarkdown

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

.PHONY: abyss abyss2 bcalm discovardenovo sga soapdenovo

abyss: \
	abyss/hsapiens-scaffolds.stats.tsv \
	abyss/GRCh38_hsapiens-scaftigs.samtobreak.tsv

abyss2: \
	abyss2/hsapiens-scaffolds.stats.tsv \
	abyss2/GRCh38_hsapiens-scaftigs.samtobreak.tsv

bcalm: \
	bcalm/hsapiens-unitigs.stats.tsv \
	bcalm/GRCh38_hsapiens-unitigs.samtobreak.tsv

discovardenovo: \
	discovardenovo/hsapiens-scaffolds.stats.tsv \
	discovardenovo/GRCh38_hsapiens-scaftigs.samtobreak.tsv

sga: \
	sga/hsapiens-contigs.stats.tsv \
	sga/GRCh38_hsapiens-contigs.samtobreak.tsv

soapdenovo: \
	soapdenovo/k63/hsapiens-scaffolds.stats.tsv \
	soapdenovo/k63/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
	soapdenovo/k79/hsapiens-scaffolds.stats.tsv \
	soapdenovo/k79/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
	soapdenovo/k95/hsapiens-scaffolds.stats.tsv \
	soapdenovo/k95/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
	soapdenovo/k111/hsapiens-scaffolds.stats.tsv \
	soapdenovo/k111/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
	soapdenovo/k127/hsapiens-scaffolds.stats.tsv \
	soapdenovo/k127/GRCh38_hsapiens-scaftigs.samtobreak.tsv

abyss_ref: \
	$(ref)-k128/$(name)-1.fa

bionano: \
	discovardenovo/besst/bionano/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
	discovardenovo/links/bionano/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
	discovardenovo/bionano/GRCh38_hsapiens-scaftigs.samtobreak.tsv

rmarkdown: \
	abyss-stats.html \
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

# Download data from NCBI.

Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz:
	curl -O ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

Homo_sapiens.GRCh38.86.chr.gff3.gz:
	curl -O ftp://ftp.ensembl.org/pub/release-86/gff3/homo_sapiens/Homo_sapiens.GRCh38.86.chr.gff3.gz

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
	sed -n 's/fastq/fq/;s/_R1//p' $< >$@

# Interleave uncorrected paired-end reads and convert lower case nucleotides to upper case
%.fq.gz: %_1.fq.gz %_2.fq.gz
	seqtk mergepe $^ | seqtk seq -U | $(gzip) >$@

# Interleave BFC-corrected paired-end reads and convert lower case nucleotides to upper case
%.bfc.fq.gz: %_1.bfc.fq.gz %_2.bfc.fq.gz
	seqtk mergepe $^ | seqtk seq -U | $(gzip) >$@

%.fa.path: %.path
	sed 's/fq.gz$$/fa.gz/' $< >$@

# Convert FASTQ to FASTA and convert lower case nucleotides to upper case
%.fa.gz: %.fq.gz
	seqtk seq -AU $< | $(gzip) >$@

# Concatenate FASTQ files
%.interleave.fq.gz: %.interleave.realpath
	cat `<$<` >$@

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

# Fill in gaps using ABySS-Sealer
%/sealer/hsapiens-scaffolds.fa: %/hsapiens-scaffolds.fa
	mkdir -p $(@D)
	time /home/benv/arch/genesis/abyss-1.9.0/k256/bin/abyss-sealer \
		-v -j$t --print-flanks \
		-L500 -F1500 -B300 \
		-o $(@D)/hsapiens -t $(@D)/hsapiens.sealer.tsv -S $< \
		-k64 -k80 -k96 -k112 -k128 -k144 -k160 -k176 -k192 -k208 -k224 -k240 \
		-ihsapiens.k64.bloom -ihsapiens.k80.bloom -ihsapiens.k96.bloom -ihsapiens.k112.bloom -ihsapiens.k128.bloom -ihsapiens.k144.bloom -ihsapiens.k160.bloom -ihsapiens.k176.bloom -ihsapiens.k192.bloom -ihsapiens.k208.bloom -ihsapiens.k224.bloom -ihsapiens.k240.bloom \
		/dev/null
	mv $(@D)/hsapiens_scaffold.fa $@

# Convert scaffolds to scaftigs
%-scaftigs.fa: %-scaffolds.fa
	seqtk seq $< | tr _ '~' | $(abyss_bin)/abyss-fatoagp -f $@ >$@.agp

# Calculate assembly contiguity metrics with abyss-fac
%.stats.tsv: %.fa
	$(abyss_bin)/abyss-fac -e3088269832 -t500 $< >$@

# Calculate assembly contiguity and correctness metrics
%.samtobreak.txt: %.sam
	(echo '==> $< <=='; bin/abyss-samtobreak-G3088269832 -l500 $<) >$@

# Convert samtobreak.txt to TSV
%.samtobreak.tsv: %.samtobreak.txt
	bin/abyss-samtobreak-to-tsv $< >$@

# Assemble the reads using ABySS 1.9
abyss/hsapiens-scaffolds.fa:
	$(MAKE) -C abyss

# Assemble the reads using ABySS 2.0
abyss2/hsapiens-scaffolds.fa:
	$(MAKE) -C abyss2

# Map the reads to the assembly using abyss-map and store the CIGAR string
k144_hsapiens-scaffolds.%.sam.cigar.gz: abyss2/k144/hsapiens-scaffolds.fa %.in
	$(abyss_bin)/abyss-map -j$t -l40 -v `<$*.in` $< | cut -sf6 | $(gzip) >$@

# BCALM

bcalm/hsapiens-unitigs.fa:
	$(MAKE) -C bcalm

# DISCOVARdenovo

discovardenovo/hsapiens-scaffolds.fa:
	$(MAKE) -C discovardenovo

# SGA

sga/hsapiens-contigs.fa:
	$(MAKE) -C sga

# SOAPdenovo

soapdenovo/k95/hsapiens-scaffolds.fa:
	$(MAKE) -C soapdenovo

# BioNano Genomics

bionano_prefix=EXP_REFINEFINAL1_q_bppAdjust_cmap_hsapiens-scaffolds_fa_NGScontigs_HYBRID_SCAFFOLD

# HybridScaffold
%/bionano/hybrid_scaffolds/$(bionano_prefix)_NOT_SCAFFOLDED.fasta %/bionano/hybrid_scaffolds/$(bionano_prefix).fasta: %/hsapiens-scaffolds.fa
	time /usr/bin/perl /gsc/btl/linuxbrew/Cellar/iryssolve/2.1.5063/scripts/HybridScaffold/hybridScaffold.pl \
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
	time bwa mem -xintractg -t$t $(ref_fa) $< >$@

# Samtools

# Sort a SAM file and output a BAM file
%.sort.bam: %.sam
	samtools sort -@$t -o $@ $<

# Index a BAM file
%.sort.bam.bai: %.sort.bam
	samtools index $<

# Calculate the genome coverage of a sorted BAM file
%.sort.bam.coverage.tsv: %.sort.bam
	samtools view -u -F260 $< | samtools depth - \
		| mlr --tsvlite --implicit-csv-header stats1 -a sum,count -f 3 \
		| mlr --tsvlite label Aligned,Covered \
		| mlr --tsvlite put '$$GenomeSize = $(GwithoutN); $$GenomeCoverage = $$Covered / $$GenomeSize' >$@

# Calculate the number of mismatches in a SAM file
%.sam.nm.tsv: %.sam
	samtools view -F260 $< | sed '/NM:i:/!d;s/^.*NM:i://;s/[[:space:]].*//' \
		| mlr --tsvlite --implicit-csv-header stats1 -a sum -f 1 \
		| mlr --tsvlite label NM >$@

# Summarize correctness and completeness
%.metrics.tsv: %.sort.bam.coverage.tsv %.sam.nm.tsv
	paste $^ | mlr --tsvlite put '$$Identity = 1 - $$NM / $$Aligned; $$QV = -10 * log10(1 - $$Identity); $$File = "$*.sam"' >$@

# QUAST

# Analayze the assemblies using QUAST
%.quast.tsv: %.fa
	quast.py -t$t -s -e -o $*.quast -R $(ref_fa) -G $(ref_gff) -l `echo $< | tr / _` $<
	cp -a $*.quast/transposed_report.tsv $@

# Assembly stats

assembly-stats.tsv: \
		abyss/k144/hsapiens-scaftigs.stats.tsv \
		abyss/k144/hsapiens-scaffolds.stats.tsv \
		abyss/k144/sealer/hsapiens-scaftigs.stats.tsv \
		abyss/k144/sealer/hsapiens-scaffolds.stats.tsv \
		abyss2/k144/hsapiens-scaftigs.stats.tsv \
		abyss2/k144/hsapiens-scaffolds.stats.tsv \
		abyss2/k144/sealer/hsapiens-scaftigs.stats.tsv \
		abyss2/k144/sealer/hsapiens-scaffolds.stats.tsv \
		bcalm/hsapiens-unitigs.stats.tsv \
		discovardenovo/hsapiens-scaftigs.stats.tsv \
		discovardenovo/hsapiens-scaffolds.stats.tsv \
		sga/hsapiens-scaftigs.stats.tsv \
		sga/hsapiens-scaffolds.stats.tsv \
		soapdenovo/k63/hsapiens-scaftigs.stats.tsv \
		soapdenovo/k63/hsapiens-scaffolds.stats.tsv \
		soapdenovo/k79/hsapiens-scaftigs.stats.tsv \
		soapdenovo/k79/hsapiens-scaffolds.stats.tsv \
		soapdenovo/k95/hsapiens-scaftigs.stats.tsv \
		soapdenovo/k95/hsapiens-scaffolds.stats.tsv \
		soapdenovo/k111/hsapiens-scaftigs.stats.tsv \
		soapdenovo/k111/hsapiens-scaffolds.stats.tsv \
		soapdenovo/k127/hsapiens-scaftigs.stats.tsv \
		soapdenovo/k127/hsapiens-scaffolds.stats.tsv \
		discovardenovo/abyss-scaffold/hsapiens-scaftigs.stats.tsv \
		discovardenovo/abyss-scaffold/hsapiens-scaffolds.stats.tsv \
		discovardenovo/besst/hsapiens-scaftigs.stats.tsv \
		discovardenovo/besst/hsapiens-scaffolds.stats.tsv \
		discovardenovo/links/hsapiens-scaftigs.stats.tsv \
		discovardenovo/links/hsapiens-scaffolds.stats.tsv \
		bionano/aggressive-B2-N2/hsapiens-scaftigs.stats.tsv \
		bionano/aggressive-B2-N2/hsapiens-scaffolds.stats.tsv \
		abyss/k144/bionano/hsapiens-scaftigs.stats.tsv \
		abyss/k144/bionano/hsapiens-scaffolds.stats.tsv \
		abyss/k144/sealer/bionano/hsapiens-scaftigs.stats.tsv \
		abyss/k144/sealer/bionano/hsapiens-scaffolds.stats.tsv \
		abyss2/k144/bionano/hsapiens-scaftigs.stats.tsv \
		abyss2/k144/bionano/hsapiens-scaffolds.stats.tsv \
		abyss2/k144/sealer/bionano/hsapiens-scaftigs.stats.tsv \
		abyss2/k144/sealer/bionano/hsapiens-scaffolds.stats.tsv \
		discovardenovo/bionano/hsapiens-scaffolds.stats.tsv \
		discovardenovo/bionano/hsapiens-scaftigs.stats.tsv \
		discovardenovo/abyss-scaffold/bionano/hsapiens-scaffolds.stats.tsv \
		discovardenovo/abyss-scaffold/bionano/hsapiens-scaftigs.stats.tsv \
		discovardenovo/besst/bionano/hsapiens-scaffolds.stats.tsv \
		discovardenovo/besst/bionano/hsapiens-scaftigs.stats.tsv \
		discovardenovo/links/bionano/hsapiens-scaffolds.stats.tsv \
		discovardenovo/links/bionano/hsapiens-scaftigs.stats.tsv
	mlr --tsvlite cat $^ >$@

samtobreak.tsv: \
		abyss/k144/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		abyss/k144/sealer/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		abyss2/k144/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		abyss2/k144/sealer/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		bcalm/GRCh38_hsapiens-unitigs.samtobreak.tsv \
		discovardenovo/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		sga/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		soapdenovo/k63/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		soapdenovo/k79/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		soapdenovo/k95/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		soapdenovo/k111/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		soapdenovo/k127/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		bionano/aggressive-B2-N2/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		abyss/k144/bionano/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		abyss/k144/sealer/bionano/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		abyss2/k144/bionano/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		abyss2/k144/sealer/bionano/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		discovardenovo/abyss-scaffold/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		discovardenovo/besst/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		discovardenovo/links/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		discovardenovo/bionano/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		discovardenovo/abyss-scaffold/bionano/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		discovardenovo/besst/bionano/GRCh38_hsapiens-scaftigs.samtobreak.tsv \
		discovardenovo/links/bionano/GRCh38_hsapiens-scaftigs.samtobreak.tsv
	mlr --tsvlite cat $^ >$@

metrics.tsv: \
		abyss/k144/GRCh38_hsapiens-scaftigs.metrics.tsv \
		abyss/k144/sealer/GRCh38_hsapiens-scaftigs.metrics.tsv \
		abyss2/k144/GRCh38_hsapiens-scaftigs.metrics.tsv \
		abyss2/k144/sealer/GRCh38_hsapiens-scaftigs.metrics.tsv \
		bcalm/GRCh38_hsapiens-unitigs.metrics.tsv \
		discovardenovo/GRCh38_hsapiens-scaftigs.metrics.tsv \
		megahit/GRCh38_hsapiens-scaftigs.metrics.tsv \
		minia/GRCh38_hsapiens-scaftigs.metrics.tsv \
		sga/GRCh38_hsapiens-scaftigs.metrics.tsv \
		soapdenovo/k95/GRCh38_hsapiens-scaftigs.metrics.tsv
	mlr --tsvlite cat $^ >$@

quast.tsv: \
		abyss/k144/hsapiens-scaffolds.quast.tsv \
		abyss2/k144/hsapiens-scaffolds.quast.tsv \
		bcalm/hsapiens-unitigs.quast.tsv \
		discovardenovo/hsapiens-scaffolds.quast.tsv \
		megahit/hsapiens-scaffolds.quast.tsv \
		minia/hsapiens-scaffolds.quast.tsv \
		sga/hsapiens-scaffolds.quast.tsv \
		soapdenovo/k95/hsapiens-scaffolds.quast.tsv
	mlr --tsvlite cat $^ >$@

# Contig stats

assembly-stats.contigs.tsv: \
		abyss/k144/hsapiens-contigs.stats.tsv \
		abyss2/k144/hsapiens-contigs.stats.tsv \
		discovardenovo/hsapiens-contigs.stats.tsv \
		sga/hsapiens-contigs.stats.tsv \
		soapdenovo/k95/hsapiens-contigs.stats.tsv
	mlr --tsvlite cat $^ >$@

samtobreak.contigs.tsv: \
		abyss/k144/GRCh38_hsapiens-contigs.samtobreak.tsv \
		abyss2/k144/GRCh38_hsapiens-contigs.samtobreak.tsv \
		discovardenovo/GRCh38_hsapiens-contigs.samtobreak.tsv \
		sga/GRCh38_hsapiens-contigs.samtobreak.tsv \
		soapdenovo/k95/GRCh38_hsapiens-contigs.samtobreak.tsv
	mlr --tsvlite cat $^ >$@

%.tsv.md: %.tsv
	mlr --itsvlite --omd cat $< >$@

# R Markdown

# Render HTML from RMarkdown
%.html: %.rmd %.tsv
	Rscript -e 'rmarkdown::render("$<", "html_document")'

assembly-stats.html: samtobreak.tsv
