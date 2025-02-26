#!/usr/bin/env bash

INPUTDIR="/data/lmattei/galapagos/24_isolates/analysis/bacass_dragonflye/bacass_out/dnaapler_output"
SCRIPTDIR="/home/ruicatxiao/fcs"
OUTPUTDIR="/data/ruicatxiao/gera_analysis/fcs_adapter_cleaned"


# Loop through all input files
find "$INPUTDIR" -name '*_reoriented.fasta' -print0 | while IFS= read -r -d '' fainput; do
    # Extract sample name from filename by removing the _reoriented.fasta suffix
    faname=$(basename "$fainput" _reoriented.fasta)
    echo "Processing sample: $faname"
    mkdir -p "${OUTPUTDIR}/${faname}"
    "$SCRIPTDIR/run_fcsadaptor.sh" \
        --fasta-input "$fainput" \
        --output-dir "${OUTPUTDIR}/${faname}" --prok \
        --container-engine singularity \
        --image "${SCRIPTDIR}/fcs-adaptor.sif"

    echo
done