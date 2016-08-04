# Assemble GIAB data using SOAPdenovo

# Parameters
k = 95

# Number of threads
t = 64

# Report run time
export SHELL=zsh -opipefail
export REPORTTIME=1
export TIMEFMT=%J  %U user %S system %P cpu %*E total %M MB

.DELETE_ON_ERROR:
.SECONDARY:

all: k$k/hsapiens-scaffolds.fa

k$k/%.scafSeq: %.config
	mkdir -p $(@D)
	cd $(@D) && time SOAPdenovo-127mer all -s ../$< -K $k -R -p $t -o $*

%-scaffolds.fa: %.scafSeq
	seqtk seq $< | tr _ '~' >$@
