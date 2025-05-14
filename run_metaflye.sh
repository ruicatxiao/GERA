#!/usr/bin/env bash

while read -r sampleid; 
do
  echo "$sampleid"
    flye \
    --nano-raw "${sampleid}_calls_filtlong.fastq.gz" \
    --out-dir "mags/${sampleid}" \
    --meta --threads 60

done < sample_id.txt

