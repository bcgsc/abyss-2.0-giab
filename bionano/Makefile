# Scaffold using BioNano Hybrid Scaffold

# Number of threads
t=32

# Reference genome
ref=GRCh38
ref_fa=/projects/btl/reference_genomes/H_sapiens/GRCh38/GCA_000001405.15_GRCh38_genomic.chr-only.fa

prefix=/gsc/btl/linuxbrew/Cellar/iryssolve/2.1.5063
bin=$(prefix)/bin
scripts=$(prefix)/scripts

scaffolds_prefix=EXP_REFINEFINAL1_q_bppAdjust_cmap_hsapiens-8_fa_NGScontigs_HYBRID_SCAFFOLD
scaffolded=$(scaffolds_prefix).fasta
not_scaffolded=$(scaffolds_prefix)_NOT_SCAFFOLDED.fasta

all: \
	default-B1-N1/hsapiens-scaftigs.samtobreak.tsv \
	default-B1-N2/hsapiens-scaftigs.samtobreak.tsv \
	default-B1-N3/hsapiens-scaftigs.samtobreak.tsv \
	default-B2-N1/hsapiens-scaftigs.samtobreak.tsv \
	default-B2-N2/hsapiens-scaftigs.samtobreak.tsv \
	default-B2-N3/hsapiens-scaftigs.samtobreak.tsv \
	default-B3-N1/hsapiens-scaftigs.samtobreak.tsv \
	default-B3-N2/hsapiens-scaftigs.samtobreak.tsv \
	default-B3-N3/hsapiens-scaftigs.samtobreak.tsv \
	aggressive-B1-N1/hsapiens-scaftigs.samtobreak.tsv \
	aggressive-B1-N2/hsapiens-scaftigs.samtobreak.tsv \
	aggressive-B1-N3/hsapiens-scaftigs.samtobreak.tsv \
	aggressive-B2-N1/hsapiens-scaftigs.samtobreak.tsv \
	aggressive-B2-N2/hsapiens-scaftigs.samtobreak.tsv \
	aggressive-B2-N3/hsapiens-scaftigs.samtobreak.tsv \
	aggressive-B3-N1/hsapiens-scaftigs.samtobreak.tsv \
	aggressive-B3-N2/hsapiens-scaftigs.samtobreak.tsv \
	aggressive-B3-N3/hsapiens-scaftigs.samtobreak.tsv

.DELETE_ON_ERROR:
.SECONDARY:

# Download data

NA24143_91x_150kb.tar.gz:
	curl -O http://www.bnxinstall.com/publicdatasets/NA24143_91x_150kb.tar.gz

# HybridScaffold

%/hybrid_scaffolds/$(scaffolded) %/hybrid_scaffolds/$(not_scaffolded): hsapiens-8.fa
	/usr/bin/perl $(scripts)/HybridScaffold/hybridScaffold.pl \
		-n $< \
		-b EXP_REFINEFINAL1_q.cmap \
		-o $* \
		-r $(bin)/RefAligner \
		-x \
		-y \
		-m all.bnx \
		-p $(scripts) \
		-q $(scripts)/optArguments_human.xml \
		-e AJmother_autoNoise1.err \
		-c `echo $* | sed 's/-/ -/g'`

%/hsapiens-scaffolds.fa: %/hybrid_scaffolds/$(not_scaffolded) %/hybrid_scaffolds/$(scaffolded)
	cat $^ >$@

# BWA

# Align scaffolds to the reference using BWA-MEM
%.sam: %.fa
	bwa mem -xintractg -t$t $(ref_fa) $< >$@

# ABySS

# Convert scaffolds to scaftigs
%-scaftigs.fa: %-scaffolds.fa
	abyss-fatoagp -f $@ $< >$@.agp

# abyss-samtobreak
%.samtobreak.txt: %.sam
	bin/abyss-samtobreak-G3088269832 -l500 $< >$@

%.samtobreak.tsv: %.samtobreak.txt
	bin/abyss-samtobreak-to-tsv $< >$@
