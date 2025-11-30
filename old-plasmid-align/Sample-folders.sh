#!/bin/bash

mkdir "INDEX"

mv *_I1_* "INDEX/"
mv *_I2_* "INDEX/"

gunzip *.fastq.gz

function folders() {
  for F in *.fastq; do
    mkdir -p Sample_`expr match "$F" '.*_S\([0-9]\{1,3\}\)'`
    mv $F Sample_`expr match "$F" '.*_S\([0-9]\{1,3\}\)'`
  done
}

folders
