# Reference genome
ref=GRCh38
ref_fa=/projects/btl/reference_genomes/H_sapiens/GRCh38/GCA_000001405.15_GRCh38_genomic.chr-only.fa

# Number of threads
t=64

# Path to ABySS executables
abyss_bin=/home/benv/arch/genesis/abyss-1.9.0/k256/bin

# Report run time and memory usage
export SHELL=zsh -opipefail
export REPORTTIME=1
export TIMEFMT=%J  %U user %S system %P cpu %*E total %M MB

.DELETE_ON_ERROR:
.SECONDARY:

all: \
	k144/$(ref)_hsapiens-scaftigs.samtobreak.tsv \
	k144/$(ref)_hsapiens-sealed-scaftigs.samtobreak.tsv \
	hsapiens-scaffolds.fa

# BWA

# Align scaffolds to the reference using BWA-MEM
$(ref)_%.sam: %.fa
	bwa mem -xintractg -t$t $(ref_fa) $< >$@

# samtools

# Sort a SAM file and output a BAM file
%.sort.bam: %.sam
	samtools sort -@$t -o $@ $<

# ABySS

k%/hsapiens-scaffolds.fa: pe400.in mp6k+unknown.in
	name=hsapiens np=$t k=$* q=15 v=-v \
		l=40 s=1000 n=10 \
		S=1000-10000 N=15 mp6k_de=--mean mp6k_n=1 \
		lib=pe400 pe400=`<pe400.in` \
		mp=mp6k mp6k=`<mp6k+unknown.in` \
		JOB_ID=k$* stages="unitigs contigs mp-sam scaffolds scaftigs" \
		TMPDIR=/var/tmp/${USER} \
		/projects/btl/bloom-dbg/assembler-comparison/common/run-abyss-1.9.0

# Construct a Bloom filter
hsapiens.k%.bloom: $(shell cat pe400.in)
	$(abyss_bin)/abyss-bloom build -v -k$* -j$t -b10G -l2 $@ $^

# Fill in gaps using ABySS-sealer
%-sealed-scaffolds.fa: %-scaffolds.fa
	$(abyss_bin)/abyss-sealer -v -j$t --print-flanks \
		-L500 -F1500 -B300 \
		-o $* -t $*.sealer.tsv -S $< \
		-k64 -k80 -k96 -k112 -k128 -k144 -k160 -k176 -k192 -k208 -k224 -k240 \
		-i$(*F).k64.bloom -i$(*F).k80.bloom -i$(*F).k96.bloom -i$(*F).k112.bloom -i$(*F).k128.bloom -i$(*F).k144.bloom -i$(*F).k160.bloom -i$(*F).k176.bloom -i$(*F).k192.bloom -i$(*F).k208.bloom -i$(*F).k224.bloom -i$(*F).k240.bloom \
		/dev/null
	mv $*_scaffold.fa $@

# Convert scaffolds to scaftigs
%-scaftigs.fa: %-scaffolds.fa
	abyss-fatoagp -f $@ $< >$@.agp

# Calculate assembly contiguity and correctness metrics
%.samtobreak.txt: %.sam
	(echo '==> $< <=='; abyss-samtobreak-G3088269832 -l500 $<) >$@

# Convert samtobreak.txt to TSV
%.samtobreak.tsv: %.samtobreak.txt
	abyss-samtobreak-to-tsv $< >$@

# Create a symlink
%-scaffolds.fa: k144/%-sealed-scaffolds.fa
	ln -sf $< $@
