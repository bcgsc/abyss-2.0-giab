# Number of threads
t=64

# Minimum overlap
m=125

# Path to ABySS executables
abyss_bin=/home/benv/arch/genesis/abyss-1.9.0/k256/bin

# Report run time and memory usage
export SHELL=zsh -opipefail
export REPORTTIME=0
export TIMEFMT=%J  %U user %S system %P cpu %*E total %M MB

.DELETE_ON_ERROR:
.SECONDARY:

all: hsapiens-contigs.fa hsapiens-scaffolds.fa

# SGA

pe400.fa.gz: %.fa.gz: %.in
	cat `<$<` >$@

mp6k+unknown.fq.gz: %.fq.gz: %.in
	cat `<$<` >$@

hsapiens.fa.gz: pe400.fa.gz
	ln -s $< $@

%.preprocess.fa: %.fa.gz
	time sga preprocess --pe-mode=2 $< >$@

%.bwt: %.fa
	time sga index -t $t -a ropebwt $<

%.filter.pass.bwt: %.fa %.bwt
	time sga filter -t $t $<

%.asqg.gz: %.bwt
	time sga overlap -t $t $*.fa -m $m

%-contigs.fa: %.preprocess.filter.pass.asqg.gz
	time sga assemble -o $* $<

# sga-align fails with reads that have been corrected using BFC.
# See https://github.com/jts/sga/issues/121

#%.refsort.bam: %.fa.gz hsapiens-contigs.fa
#	time sga-align -t $t --name $* hsapiens-contigs.fa $<

#%.refsort.bam: %.fq.gz hsapiens-contigs.fa
#	time sga-align -t $t --name $* hsapiens-contigs.fa $<

%.bam: hsapiens-contigs.fa %.fa.gz
	time bwa mem -t$t -p $^ | samtools view -F2304 -b -o $@ -

%.bam: hsapiens-contigs.fa %.fq.gz
	time bwa mem -t$t -p $^ | samtools view -F2304 -b -o $@ -

%.refsort.bam: %.bam
	time samtools sort -@$t -o $@ $<

%.de: %.hist %.diffcontigs.sorted.bam
	time DistanceEst -v -s 500 --mind -99 -n 10 -k 99 -j $t -l 100 -o $@ $^

%.de: %.bam
	time sga-bam2de.pl -t $t -m 500 -n 10 --prefix $* $<

%.astat: %.refsort.bam
	time sga-astat.py -m 500 $< >$@

%.scaf: %-contigs.fa %-graph.asqg.gz pe400.astat pe400.de mp6k+unknown.de
	time sga scaffold -m 500 -g $*-graph.asqg.gz -a pe400.astat --pe pe400.de --mate-pair mp6k+unknown.de -o $@ $<

# sga scaffold2fasta blows an assertion.
# See https://github.com/jts/sga/issues/125
#%-scaffolds.fa: %.scaf %-graph.asqg.gz
#	time sga scaffold2fasta -m 500 --write-unplaced --use-overlap -o $@ -a $*-graph.asqg.gz $<

%.path: %.scaf
	../bin/sga-scaf-to-path $< >$@

%-scaffolds.fa: %-contigs.fa %.path
	$(abyss_bin)/MergeContigs -v -k100 -o $@ $< $< $*.path
