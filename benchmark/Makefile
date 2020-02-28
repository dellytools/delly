SHELL := /bin/bash

# Targets
TARGETS = .conda .channels .envs .truvari .delly .visor
PBASE=$(shell pwd)

all:   	$(TARGETS)

.conda:
	wget 'https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh' && bash Miniconda3-latest-Linux-x86_64.sh -b -p ${PBASE}/bin && rm -f Miniconda3-latest-Linux-x86_64.sh && touch .conda

.channels: .conda
	export PATH=${PBASE}/bin/bin:${PATH} && conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda && touch .channels

.envs: .conda .channels
	export PATH=${PBASE}/bin/bin:${PATH} && conda create -y --prefix=${PBASE}/bin/envs/sv minimap2 samtools bcftools bedtools htslib wgsim pbsim bwa alfred && conda create -y --prefix=${PBASE}/bin/envs/truvari samtools htslib && touch .envs

.truvari: .conda .channels .envs
	export PATH=${PBASE}/bin/bin:${PATH} && source activate ${PBASE}/bin/envs/truvari && pip install biopython truvari==1.2 && conda deactivate && touch .truvari

.visor: .conda .channels .envs
	export PATH=${PBASE}/bin/bin:${PATH} && source activate ${PBASE}/bin/envs/sv && git clone https://github.com/davidebolo1993/VISOR.git && cd VISOR && python setup.py install && conda deactivate && cd ../ && touch .visor

.delly: .conda .channels .envs .truvari
	cd ../ && make all && make install && cd benchmark && touch .delly

test: .conda .channels .envs .truvari .delly
	export PATH=${PBASE}/bin/bin:${PATH} && source activate sv && VISOR -h && ../bin/dellyLR --version && svim --version && conda deactivate

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) bin/ VISOR/ sim_*/ Miniconda3-latest-Linux-x86_64.sh*
