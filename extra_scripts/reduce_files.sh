#!/bin/bash
# script reduces the fastq files to make a subset on which the pipeline can be tested within 10 minutes

for f in /nfs/data/Big_Data/L15-217_framboos/data/raw_sequences_names_2/*
do
  echo "Processing $f file..."
  # take action on each file. $f store current file name
  full_name="${f##*/}_reduced" # filename
  extention=".${f#*.}" # extention
  zcat $f | head -n 51010 | bgzip -c > "/nfs/data/Big_Data/L15-217_framboos/data/raw_sequences_reduced/$full_name$extention"
done
