# Number of threads
t=64

# Minimum overlap
m=125

# Report run time and memory usage
export SHELL=zsh -opipefail
export REPORTTIME=0
export TIMEFMT=%J  %U user %S system %P cpu %*E total %M MB

.DELETE_ON_ERROR:
.SECONDARY:

all: hsapiens-contigs.fa
#all: testreads-contigs.fa

# SGA

hsapiens.fa.gz: pe400.in
	cat `<$<` >$@

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
