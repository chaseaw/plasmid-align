#!/bin/bash

# Setup Ubuntu
sudo apt update --yes
sudo apt upgrade --yes

# check if conda installed
conDIR=`which conda`
conVAL=${#conDIR}

# if conda is not installed
if [ "$conVAL" -eq "0" ]; then

	# Get Miniconda3
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
	bash ~/miniconda.sh -b -p ~/miniconda3 
	rm ~/miniconda.sh

	# Source terminal so conda can be run
	FILE=~/.bashrc
	if test -f "$FILE"; then
    	source $FILE
	fi

	FILE=~/.bash_profile
	if test -f "$FILE"; then
    	source $FILE
	fi
else
	echo "conda installation detected=$conDIR"
fi

# update conda
yes | conda update conda 

# add expanded package directories
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge 

# create conda environment for plasmid alignment
conda create -n plasmid-align bwa-mem2 samtools megahit

exec $SHELL
