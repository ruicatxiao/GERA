#!/usr/bin/env bash

SOURCEDIR="/data/ruicatxiao/gera_analysis/amr"
# GENOMEDIR="/data/ruicatxiao/gera_analysis/genome"
PLASMIDIR="/data/ruicatxiao/gera_analysis/plasmid"
OUTDIR="/data/ruicatxiao/gera_analysis/amr"

find "$SOURCEDIR" -name '*_genome.fasta' -print0 | while IFS= read -r -d '' fainput; do
    faname=$(basename "$fainput" _genome.fasta)
    echo
    echo "========================="
    echo "${faname} GENOME IS BEING PROCESSED"

    amrfinder \
        -p "${SOURCEDIR}/${faname}_genome_prot.fasta" \
        -g "${SOURCEDIR}/${faname}_genome_amr.gff" \
        -n "${SOURCEDIR}/${faname}_genome.fasta" \
        -O Escherichia --plus \
        -o "${OUTDIR}/${faname}_genome_amrplus.tsv"

    echo "${faname} GENOME DONE"
    echo
    echo "${faname} PLASMID IS BEING PROCESSED"

    amrfinder \
        -p "${SOURCEDIR}/${faname}_plasmid_prot.fasta" \
        -g "${SOURCEDIR}/${faname}_plasmid_amr.gff" \
        -n "${PLASMIDIR}/${faname}_plasmid.fasta" \
        -O Escherichia --plus \
        -o "${OUTDIR}/${faname}_plasmid_amrplus.tsv"

    echo "${faname} PLASMID DONE"
    echo "========================="
done