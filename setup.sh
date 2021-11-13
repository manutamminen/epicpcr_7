#!/bin/bash

# Download and prepare raw data
mkdir -p data/raw
cd data/raw
cat ../../zenodo_files.txt |
  while read file;
  do wget $file;
  done
cd ../..

# Archive the old final output file
cd docs
mv index.md archived-report-md
cd ..
