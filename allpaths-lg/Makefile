all: hsapiens-scaffolds.fa

# Report run time and memory usage
export SHELL=zsh -opipefail
export REPORTTIME=1
export TIMEFMT=%J  %U user %S system %P cpu %*E total %M MB

.DELETE_ON_ERROR:
.SECONDARY:

hsapiens-scaffolds.fa:
	PrepareAllPathsInputs.pl DATA_DIR=$(PWD) PLOIDY=2 HOSTS=32
	RunAllPathsLG PRE=. REFERENCE_NAME=. DATA_SUBDIR=. RUN=allpaths SUBDIR=run
