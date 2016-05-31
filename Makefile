# Assemble Genome in a Bottle
# Written by Shaun Jackman <sjackman@gmail.com>

# Name of the species and assemblies
name=hsapiens

# Number of threads
t=8

# Reference genome
ref=GRCh38
ref_fa=/genesis/extscratch/btl/reference_genomes/H_sapiens/GRCh38/GCA_000001405.15_GRCh38_genomic.chr-only.fa
ref_gff=/genesis/extscratch/btl/reference_genomes/H_sapiens/GRCh38/Homo_sapiens.GRCh38.84.chr.gff

# Report run time
export SHELL=zsh -opipefail
export REPORTTIME=1

.DELETE_ON_ERROR:
.SECONDARY:

all: curl kmerstream nxtrim bfc bfc_kmerstream abyss_ref quast rmarkdown

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

abyss_ref: \
	$(ref)-k128/$(name)-1.fa

quast: \
	abyss/k96/$(name)-1.quast/transposed_report.tsv \
	abyss/k128/$(name)-1.quast/transposed_report.tsv \
	abyss/k160/$(name)-1.quast/transposed_report.tsv \
	abyss/k192/$(name)-1.quast/transposed_report.tsv

rmarkdown: \
	assembly-stats.html

# Install dependencies

deps=bfc fastqc jq kmerstream nxtrim

deps:
	@echo $(deps)

install-deps:
	brew install $(deps)

versions.json:
	brew info --json=v1 $(deps) >$@

versions.tsv: %.tsv: %.json
	jq -r '["Name", "Version"],(.[]|[.name,.linked_keg])|@tsv' $< >$@

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
%.kmerstream: %.realpath
	KmerStream -t$t -s2 -k31,32,33,63,64,65,95,96,97,127,128,129,159,160,161,191,192,193,223,224,225 -o$@ `<$<`

# Convert to TSV
%.kmerstream.tsv: %.kmerstream
	(printf "Q\tk\tF0\tf1\tF1\n"; \
		gsed 's/[^ ]* = //g;s/, /\n/' $< | paste -d'\t' - - - - -) \
		| estimate.py /dev/stdin >$@

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
	bfc -t$t -s3G -r $(bfc) /dev/null $< | gzip >$@

# Correct the reads
%.bfc.log: %.path $(bfc)
	$(MAKE) `sed 's/^/$(<D)\//;s/.fastq.gz/.bfc.fq.gz/' $<` 2>&1 | tee $@

%.bfc.path: %.path
	sed 's/.fastq.gz/.bfc.fq.gz/' $< >$@

HG004/mp6k.mp.path:
	printf "MPHG004-23100079/MPHG004_S3_L003_001.mp.fastq.gz\nMPHG004-23110109/MPHG004_S3_L003_001.mp.fastq.gz\n" >$@

HG004/mp6k.unknown.path:
	printf "MPHG004-23100079/MPHG004_S3_L003_001.unknown.fastq.gz\nMPHG004-23110109/MPHG004_S3_L003_001.unknown.fastq.gz\n" >$@

# ABySS

# Assemble the reference genome
$(ref)-k%/$(name)-1.fa: $(ref_fa)
	mkdir -p $(ref)-k$*
	ABYSS -v -k$* -e0 -t0 -c0 $< -o $@ -s $(ref)-k$*/$(name)-bubbles.fa 2>&1 |tee $@.log

# QUAST

# Analayze the assemblies using QUAST
%.quast/transposed_report.tsv: %.fa
	quast.py -t$t -eL -o $*.quast -R $(ref_fa) -G $(ref_gff) $<

# R Markdown

# Render HTML from RMarkdown
%.html: %.rmd %.tsv
	Rscript -e 'rmarkdown::render("$<", "html_document")'
