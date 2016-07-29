# Assemble GIAB data using SOAPdenovo

# Parameters
k = 63

# Number of threads
t = 64

# Report run time
export SHELL=zsh -opipefail
export REPORTTIME=1
export TIMEFMT=%J  %U user %S system %P cpu %*E total %M MB

.DELETE_ON_ERROR:
.SECONDARY:

all: hsapiens-scaffolds.fa

%.scafSeq: %.config
	time SOAPdenovo-63mer all -s $< -K $k -R -p $t -o $*

%-scaffolds.fa: %.scafSeq
	seqtk seq $< | tr _ '~' >$@
