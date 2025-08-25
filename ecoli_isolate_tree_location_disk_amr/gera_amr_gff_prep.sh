#!/usr/bin/env bash

SOURCEFADIR="gera_analysis/amr"
SOURCEGFFDIR="gera_analysis/genome"
OUTDIR="gera_analysis/amr"

find "$SOURCEFADIR" -name '*_genome.fasta' -print0 | while IFS= read -r -d '' fainput; do
    faname=$(basename "$fainput" _genome.fasta)

    cat "${SOURCEGFFDIR}/${faname}_genome.gff" | awk '$1!~"^#"' \
    | awk '$3=="gene" || $3=="pseudogene"' \
    | awk -F '\t' '{split($9,a,";"); print $1, ".", "gene", $4, $5, $6, $7, $8, a[1]}' OFS='\t' \
    | sed 's/ID=/Name=/g' | sed 's/agat-gene-/agat-rna-/g' > "${OUTDIR}/${faname}_genome_amr.gff"

done