#!/bin/bash

INPUT_FOLDER="/path/to/your/input"
OUTPUT_FOLDER="/path/to/your/output"

mkdir -p "${OUTPUT_FOLDER}"/{bam,demuxed,fastq}

for run_dir in $(find "$INPUT_FOLDER" -maxdepth 2 -type d -name "*_run_*"); do
    run_name=$(basename "$run_dir")
    pod5_dir="${run_dir}/pod5"
    
    [[ ! -d "$pod5_dir" ]] && continue
    
    bam_output="${OUTPUT_FOLDER}/bam/${run_name}.bam"
    
    dorado basecaller sup,4mC_5mC,6mA "$pod5_dir" > "$bam_output" --device cuda:0,1
    dorado demux --kit-name SQK-NBD114-24 --output-dir "${OUTPUT_FOLDER}/demuxed" "$bam_output"
    
    for bam_file in "${OUTPUT_FOLDER}/demuxed/${run_name}_barcode"*.bam; do
        [[ ! -f "$bam_file" ]] && continue
        bam_basename=$(basename "$bam_file" .bam)
        bedtools bamtofastq -i "$bam_file" -fq "${OUTPUT_FOLDER}/fastq/${bam_basename}.fastq"
    done
done