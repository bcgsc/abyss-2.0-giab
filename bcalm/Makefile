# Number of threads
t=64

# Size of a k-mer
k=63

# Report run time and memory usage
export SHELL=zsh -opipefail
export REPORTTIME=1
export TIMEFMT=%J  %U user %S system %P cpu %*E total %M MB

.DELETE_ON_ERROR:
.SECONDARY:

all: hsapiens-unitigs.fa

# BCALM

hsapiens-unitigs.h5: pe400.in
	time bcalm -in $< -out hsapiens-unitigs -k $k -abundance 5 -nb-cores $t

%.fa: %.h5
	time bglue -in $< -out $@ -k $k -nb-cores 12
