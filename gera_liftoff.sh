#!/usr/bin/env bash

TARGETFADIR="/data/lmattei/galapagos/24_isolates/analysis/bacass_dragonflye/bacass_out/dnaapler_output"
SOURCEFADIR="/data/lmattei/galapagos/24_isolates/analysis/bacass_dragonflye/bacass_out/Medaka"

# Not using this source gff, as I have used agat to clean it already
# SOURCEGFFDIR="/data/lmattei/galapagos/24_isolates/analysis/bacass_dragonflye/bacass_out/Prokka"

OUTPUTDIR="/data/ruicatxiao/gera_analysis/liftoff_fa_gff"


# Loop through all input files
find "$TARGETFADIR" -name '*_reoriented.fasta' -print0 | while IFS= read -r -d '' fainput; do
    # Extract sample name from filename by removing the _reoriented.fasta suffix
    echo
    echo "===START A NEW SAMPLE===" 
    faname=$(basename "$fainput" _reoriented.fasta)
    echo "Target fasta sample name: $faname"
    ls "${TARGETFADIR}/${faname}/${faname}_reoriented.fasta"
    echo
    echo "Source fasta sample name: ${SOURCEFADIR}/${faname}_polished_genome.fa/consensus.fasta"
    ls "${SOURCEFADIR}/${faname}_polished_genome.fa/consensus.fasta"
    echo
    # echo "Source gff sample name: ${SOURCEGFFDIR}/${faname}/${faname}.gff"
    # ls "${SOURCEGFFDIR}/${faname}/${faname}.gff"

    liftoff \
        "${TARGETFADIR}/${faname}/${faname}_reoriented.fasta" \
        "${SOURCEFADIR}/${faname}_polished_genome.fa/consensus.fasta" \
        -g "${OUTPUTDIR}/${faname}_cleaned.gff" \
        -o "${OUTPUTDIR}/${faname}.gff"

    cp "${TARGETFADIR}/${faname}/${faname}_reoriented.fasta" "${OUTPUTDIR}/${faname}.fasta"

    echo
    echo "===END A SAMPLE==="
    echo

done