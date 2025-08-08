#!/bin/bash

mkdir -p dnaapler_output

for file in *.fasta; do
    sample_name=$(basename "$file" .fasta)
    dnaapler all -i "$file" -o "dnaapler_output/${sample_name}" -f -p "$sample_name" -t 10
    if [ $? -eq 0 ]; then
        echo "Processed: $file -> dnaapler_output/${sample_name}"
    else
        echo "Failed to process: $file"
    fi
done