SHELL := /bin/bash

# Targets
TARGETS = .conda .channels .envs .truvari .delly
PBASE=$(shell pwd)

all:   	$(TARGETS)

.conda:
	wget 'https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh' && bash Miniconda3-latest-Linux-x86_64.sh -b -p ${PBASE}/bin && rm -f Miniconda3-latest-Linux-x86_64.sh && touch .conda

.channels: .conda
	export PATH=${PBASE}/bin/bin:${PATH} && conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda && touch .channels

.envs: .conda .channels
	export PATH=${PBASE}/bin/bin:${PATH} && conda create -y --prefix=${PBASE}/bin/envs/sv svim sniffles samtools bcftools htslib && touch .envs

.truvari: .conda .channels .envs
	export PATH=${PBASE}/bin/bin:${PATH} && source activate ${PBASE}/bin/envs/sv && pip install truvari && conda deactivate && touch .truvari

.delly: .conda .channels .envs .truvari
	cd ../ && make all && make install && cd benchmark && touch .delly

test: .conda .channels .envs .truvari .delly
	export PATH=${PBASE}/bin/bin:${PATH} && source activate sv && ../bin/dellyLR --version && svim --version && sniffles --help && truvari --help && source deactivate

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) bin/
