#!/bin/bash

function folderz() {
  for F in *.fastq.gz; do
    directory="Sample_$(echo $F | grep -oP '.*_S\K[0-9]{1,3}')"
    file_prefix="${F%_S[0-9]*_R[0-9]*_001.fastq.gz}"
    mkdir -p "$directory"
    mv "$F" "$directory"
    echo "$directory/$file_prefix" >> fastq_folders.csv
  done
}

folderz

sort -u fastq_folders.csv -o fastq_folders.csv
