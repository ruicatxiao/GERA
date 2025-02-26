#!/usr/bin/env bash

SOURCEFADIR="/data/ruicatxiao/gera_analysis/amr"
SOURCEGFFDIR="/data/ruicatxiao/gera_analysis/genome"
OUTDIR="/data/ruicatxiao/gera_analysis/amr"




# Loop through all input files, this part for genome protein
find "$SOURCEFADIR" -name '*_genome.fasta' -print0 | while IFS= read -r -d '' fainput; do
    echo
    echo "===START A NEW SAMPLE===" 
    faname=$(basename "$fainput" _genome.fasta)

    docker run --rm \
            -v "$SOURCEFADIR:/input_fa" \
            -v "$SOURCEGFFDIR:/input_gff" \
            -v "$OUTDIR:/output" \
            quay.io/biocontainers/agat:1.4.0--pl5321hdfd78af_0 \
            agat_sp_extract_sequences.pl \
            -f "/input_fa/${faname}_genome.fasta" \
            -g "/input_gff/${faname}_genome.gff" \
            -t cds -p \
            --output "/output/${faname}_genome_prot.fasta"

    echo
    echo "===END A SAMPLE==="
    echo

done