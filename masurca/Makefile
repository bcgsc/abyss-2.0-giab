all: hsapiens-scaffolds.fa

.DELETE_ON_ERROR:
.SECONDARY:

assemble.sh: config.txt
	masurca $<

hsapiens-scaffolds.fa: assemble.sh
	./assemble.sh
