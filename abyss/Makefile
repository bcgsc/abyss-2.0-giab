# Reference genome
ref=GRCh38
ref_fa=/genesis/extscratch/btl/reference_genomes/H_sapiens/GRCh38/GCA_000001405.15_GRCh38_genomic.chr-only.fa

# Number of threads
t=12

# Path to ABySS executables
abyss_bin=/home/benv/arch/genesis/abyss-1.9.0/k256/bin

.DELETE_ON_ERROR:
.SECONDARY:

all: \
	k96/$(ref)_hsapiens-sealed-scaftigs.sort.bam \
	k112/$(ref)_hsapiens-sealed-scaftigs.sort.bam \
	k128/$(ref)_hsapiens-sealed-scaftigs.sort.bam \
	k144/$(ref)_hsapiens-sealed-scaftigs.sort.bam \
	k160/$(ref)_hsapiens-sealed-scaftigs.sort.bam \
	k176/$(ref)_hsapiens-sealed-scaftigs.sort.bam \
	k192/$(ref)_hsapiens-sealed-scaftigs.sort.bam

bwamem: \
	k96/$(ref)_hsapiens-scaftigs.sort.bam \
	k112/$(ref)_hsapiens-scaftigs.sort.bam \
	k128/$(ref)_hsapiens-scaftigs.sort.bam \
	k136/$(ref)_hsapiens-scaftigs.sort.bam \
	k140/$(ref)_hsapiens-scaftigs.sort.bam \
	k144/$(ref)_hsapiens-scaftigs.sort.bam \
	k148/$(ref)_hsapiens-scaftigs.sort.bam \
	k152/$(ref)_hsapiens-scaftigs.sort.bam \
	k160/$(ref)_hsapiens-scaftigs.sort.bam \
	k176/$(ref)_hsapiens-scaftigs.sort.bam \
	k192/$(ref)_hsapiens-scaftigs.sort.bam

# Align scaffolds to the reference using BWA-MEM
$(ref)_%.sam: %.fa
	bwa mem -xintractg -t$t $(ref_fa) $< >$@

# Sort a SAM file and output a BAM file
%.sort.bam: %.sam
	samtools sort -@$t -o $@ $<

# Construct a Bloom filter
hsapiens.k%.bloom: $(shell cat pe400.in)
	$(abyss_bin)/abyss-bloom build -v -k$* -j$t -b10G -l2 $@ $^

# Fill in gaps using ABySS-sealer
%-sealed-scaffolds.fa: %-scaffolds.fa
	$(abyss_bin)/abyss-sealer -v -j$t --print-flanks \
		-L500 -F1500 -B300 \
		-o $* -t $*.sealer.tsv -S $< \
		-k64 -k96 -k128 -k160 -k192 -k224 \
		-i$(*F).k64.bloom -i$(*F).k96.bloom -i$(*F).k128.bloom -i$(*F).k160.bloom -i$(*F).k192.bloom -i$(*F).k224.bloom \
		/dev/null
	mv $*_scaffold.fa $@

# Convert scaffolds to scaftigs
%-scaftigs.fa: %-scaffolds.fa
	abyss-fatoagp -f $@ $< >$@.agp
