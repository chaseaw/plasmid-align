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
	FILE=`echo ~/.bashrc`
	if [[ -f "$FILE" ]]; then
		echo 'export PATH="~/miniconda3/bin:$PATH"' >> $FILE
    	source $FILE
	fi

	FILE=`echo ~/.bash_profile`
	if [[ -f "$FILE" ]]; then
		echo 'export PATH="~/miniconda3/bin:$PATH"' >> $FILE
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
conda create -n plasmid-align bwa-mem2=2.1 samtools=1.11 megahit=1.2.9 trimmomatic=0.39

# initialize conda
conda init bash

exec $SHELL

