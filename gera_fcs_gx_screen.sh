#!/usr/bin/env bash

INPUTDIR="/data/ruicatxiao/gera_analysis/fcs_adapter_cleaned/cleaned_sequences"
SCRIPTDIR="/home/ruicatxiao/fcs"
OUTPUTDIR="/data/ruicatxiao/gera_analysis/fcs_gx_cleaned"
GXDBDIR="/home/ruicatxiao/fcs/fcs_gxdb"
TAXID=562

# Loop through all input files
find "$INPUTDIR" -name '*_reoriented.fasta' -print0 | while IFS= read -r -d '' fainput; do
    # Extract sample name from filename by removing the _reoriented.fasta suffix
    faname=$(basename "$fainput" _reoriented.fasta)
    echo "Processing sample: $faname"
    python3 "$SCRIPTDIR/fcs.py" \
        screen genome \
        --fasta "$fainput" \
        --out-dir "$OUTPUTDIR" \
        --gx-db "$GXDBDIR" --tax-id "$TAXID" < /dev/null
    echo
done