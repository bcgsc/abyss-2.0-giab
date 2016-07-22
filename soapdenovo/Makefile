# Assemble GIAB data using SOAPdenovo

# Parameters
k = 63

# Number of threads
t = 32

# Genome size
G = 3088269832

# Report run time
export SHELL=zsh -opipefail
export REPORTTIME=1

.DELETE_ON_ERROR:
.SECONDARY:

all: hsapiens.scafSeq

%.scafSeq: %.config
	SOAPdenovo-63mer all -s $< -K$k -R -p$t -N$G -o $*
	ln -sf $@ $*-scaffolds.fa
