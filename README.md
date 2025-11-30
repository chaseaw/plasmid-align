Reference/De Novo Alignment of Illumina/Nanopore reads
================================================================================
Description: This directory will enable easy analysis of sequenced plasmids and amplicons.

There are a number of ways to do this but this one should be comprehensive and relatively painless.

[Go To Installation](#installation)

[Go To Usage](#usage)

Prerequisites for using this pipeline
================================================================================

**1.** **Git**

*  From your terminal:

```
sudo apt update
sudo apt install git
```

**2.** **IGV** (for viewing alignments)

Get IGV for your operating system (you don't need the command line version):
* [IGV Software Download](http://software.broadinstitute.org/software/igv/download)

**3.** **Conda and Dependencies** (see conda yml for details)
* python-based scipt
* fastp for trimming
* bwa aligner for illumina reference alignment
* minimap2 for nanopore reference alignment
* spades for illumina de novo alignment
* flye for nanopore de novo alignment
* samtools and seqtk

Installation
================================================================================
From your Software directory on the command line:

```
git clone https://github.com/chaseaw/plasmid_seq.git
```

This will create a directory called plasmid_seq that includes all necessary scripts.

To create replica conda environment on great lakes cluster:
```
conda create --name plasmid_seq --file plasmid_seq-spec.txt
```
To assemble replica environment on other systems:
```
conda env create -f plasmid_seq.yml -n plasmid_seq
```
To assemble environment from scratch:
```
conda create -n plasmid_seq -c conda-forge -c bioconda -c defaults \
    python=3.11 \
    bwa \
    fastp \
    flye \
    minimap2 \
    samtools \
    spades \
    seqtk
```
Add plasmid_seq.py and make executable to your path for future use.

Usage
================================================================================

```
Example usages:
  # Illumina reads aligned to a plasmid or linear amplicon reference
  python plasmid_seq.py --illumina sample_R1.fq.gz sample_R2.fq.gz \
              --reference plasmid.fasta \
              --output plasmid_align

  # Nanopore reads aligned to a reference (plasmid or amplicon)
  python plasmid_seq.py --nanopore nanopore.fastq \
              --reference plasmid.fasta \
              --output ont_align

  # Illumina de novo assembly
  python plasmid_seq.py --illumina R1.fq.gz R2.fq.gz \
              --de_novo \
              --output illumina_assembly

  # Nanopore de novo assembly
  python plasmid_seq.py --nanopore longreads.fq \
              --de_novo \
              --output ont_assembly
```

Additional flags:

--threads specifies how many proccessors to dedicate to the run

--genome_size is explicitly for de novo plasmid alignment with Nanopore reads (default is "15k")

Notes:

De Novo alignment of nanopore reads will need more than 8GB RAM to not throw an error. 

Run IGV loading the reference .fa as genome and your new .bam as a track to identify any sequence changes.
