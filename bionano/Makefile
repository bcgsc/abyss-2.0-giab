# Scaffold using BioNano Hybrid Scaffold

# Number of threads
t=32

# Report run time and memory usage
export SHELL=zsh -opipefail
export REPORTTIME=1
export TIMEFMT=%J  %U user %S system %P cpu %*E total %M MB

# Path to ABySS executables
abyss_bin=/gsc/btl/linuxbrew/Cellar/abyss/1.9.0-k256/bin

# Reference genome
ref=GRCh38
ref_fa=/projects/btl/reference_genomes/H_sapiens/GRCh38/GCA_000001405.15_GRCh38_genomic.chr-only.fa

prefix=/gsc/btl/linuxbrew/Cellar/iryssolve/2.1.5063
bin=$(prefix)/bin
scripts=$(prefix)/scripts

scaffolds_prefix=EXP_REFINEFINAL1_q_bppAdjust_cmap_hsapiens-8_fa_NGScontigs_HYBRID_SCAFFOLD
scaffolded=$(scaffolds_prefix).fasta
not_scaffolded=$(scaffolds_prefix)_NOT_SCAFFOLDED.fasta

all: aggressive-B2-N2/hsapiens-sealed-scaftigs.samtobreak.tsv

optimize-hybridscaffold: \
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

# Fill in gaps using ABySS-sealer
%-sealed-scaffolds.fa: %-scaffolds.fa
	$(abyss_bin)/abyss-sealer -v -j$t --print-flanks \
		-L500 -F1500 -B300 \
		-o $* -t $*.sealer.tsv -S $< \
		-k64 -k80 -k96 -k112 -k128 -k144 -k160 -k176 -k192 -k208 -k224 -k240 \
		-i$(*F).k64.bloom -i$(*F).k80.bloom -i$(*F).k96.bloom -i$(*F).k112.bloom -i$(*F).k128.bloom -i$(*F).k144.bloom -i$(*F).k160.bloom -i$(*F).k176.bloom -i$(*F).k192.bloom -i$(*F).k208.bloom -i$(*F).k224.bloom -i$(*F).k240.bloom \
		/dev/null
	mv $*_scaffold.fa $@

# Fill in gaps using ABySS-sealer varying parameter F
dir=aggressive-B2-N2
name=hsapiens
$(dir)/F%/$(name)-sealer-scaffolds.fa: $(dir)/$(name)-scaffolds.fa
	mkdir $(@D)
	$(abyss_bin)/abyss-sealer -v -j$t --print-flanks \
		-L500 -F$* -B300 \
		-o $(@D)/$(name) -t $(@D)/$(name).sealer.tsv -S $< \
		-k64 -k80 -k96 -k112 -k128 -k144 -k160 -k176 -k192 -k208 -k224 -k240 \
		-i$(name).k64.bloom -i$(name).k80.bloom -i$(name).k96.bloom -i$(name).k112.bloom -i$(name).k128.bloom -i$(name).k144.bloom -i$(name).k160.bloom -i$(name).k176.bloom -i$(name).k192.bloom -i$(name).k208.bloom -i$(name).k224.bloom -i$(name).k240.bloom \
		/dev/null
	mv $(@D)/$(name)_scaffold.fa $@

# Convert scaffolds to scaftigs
%-scaftigs.fa: %-scaffolds.fa
	abyss-fatoagp -f $@ $< >$@.agp

# abyss-samtobreak
%.samtobreak.txt: %.sam
	(echo '==> $< <=='; bin/abyss-samtobreak-G3088269832 -l500 $<) >$@

%.samtobreak.tsv: %.samtobreak.txt
	bin/abyss-samtobreak-to-tsv $< >$@
