#!/usr/bin/env bash

TARGETFADIR="bacass_dragonflye/bacass_out/dnaapler_output"
SOURCEFADIR="bacass_dragonflye/bacass_out/Medaka"
SOURCEGFFDIR="bacass_dragonflye/bacass_out/Prokka"
OUTPUTDIR="gera_analysis/liftoff_fa_gff"


# Loop through all input files
find "$TARGETFADIR" -name '*_reoriented.fasta' -print0 | while IFS= read -r -d '' fainput; do
    # Extract sample name from filename by removing the _reoriented.fasta suffix
    echo
    echo "===START A NEW SAMPLE===" 
    faname=$(basename "$fainput" _reoriented.fasta)

    echo "Source gff sample name: ${SOURCEGFFDIR}/${faname}/${faname}.gff"
    ls "${SOURCEGFFDIR}/${faname}/${faname}.gff"

    docker run --rm \
            -v "$SOURCEGFFDIR:/input_gff" \
            -v "$OUTPUTDIR:/output" \
            quay.io/biocontainers/agat:1.4.0--pl5321hdfd78af_0 \
            agat_convert_sp_gxf2gxf.pl \
            -g "/input_gff/${faname}/${faname}.gff" \
            -o "/output/${faname}_cleaned.gff"

    echo
    echo "===END A SAMPLE==="
    echo

done