#!/bin/bash

#AMR Finder Plus on final nanopore metagenome assemblies

INPUT_DIR="/data/lmattei/galapagos/anvio/02_assembly"   \
OUTPUT_DIR="/data/lmattei/galapagos/amrfinderplus/out"

mkdir -p "$OUTPUT_DIR"

# Loop through .fasta files in INPUT_DIR
for assembly in "$INPUT_DIR"/*.fasta; do

    # Extract the sample name (e.g., "22_baquerizo") from the filename
    sample_name=$(basename "$assembly" .fasta | cut -d'_' -f3-)

    # Define output file path
    output_file="$OUTPUT_DIR/${sample_name}.amrfinder"

    # Run AMRFinder+
    echo "Processing: $(basename "$assembly") -> $(basename "$output_file")"
    amrfinder -n "$assembly" --threads 8 --plus --name="consensus_proovframe_${sample_name}" \
        > "$output_file"
done

head -1 $(ls *.amrfinder | head -1) > amrfinder_combined.tsv
grep -h -v 'Name' *.amrfinder >> amrfinder_combined.tsv

#manually remove stress and virulence genes > amrfinder_combined_cleaned.tsv
