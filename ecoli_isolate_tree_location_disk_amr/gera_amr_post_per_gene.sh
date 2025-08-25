#!/usr/bin/env bash
set -euo pipefail

INPUTDIR="gera_analysis/amr"
SAMPLE_LIST="$INPUTDIR/amr_sampleid.txt"
OUTPUT="$INPUTDIR/amr_combined_summary.tsv"
# GENE SYMBOLS
gene_symbols=(
  "aph(3')-Ia" "blaTEM-176" "blaTEM" "blaTEM-1" "emrD" "floR" "acrF"
  "blaCTX-M-3" "blaCTX-M-15" "blaCTX-M-32" "blaCTX-M-65" "blaEC" "erm(B)"
  "mph(A)" "gyrA_D87N" "gyrA_S83A" "gyrA_S83L" "parC_S57T" "parC_S80I"
  "parE_S458A" "qnrB19" "qnrS1" "mdtM" "aadA1" "aadA5" "aph(3'')-Ib"
  "aph(6)-Id" "tet(A)" "tet(B)" "dfrA1" "dfrA7" "dfrA14" "dfrA17"
  "dfrA51" "sul1" "sul2" "sul3" "cyaA_S352T" "fosA3" "sat2" "uhpT_E350Q"
)
{
  printf "SAMPLEID"
  for sym in "${gene_symbols[@]}"; do
    printf "\t%s" "$sym"
  done
  printf "\n"
} > "$OUTPUT"

while IFS= read -r sample || [ -n "$sample" ]; do
  sample=$(echo "$sample" | xargs)
  GENOME_FILE="$INPUTDIR/${sample}_genome_amrplus.tsv"
  PLASMID_FILE="$INPUTDIR/${sample}_plasmid_amrplus.tsv"
  declare -A gset=() pset=()
  if [[ -f "$GENOME_FILE" ]]; then
    while IFS=$'\t' read -r _ contig start stop strand symbol _ _ type _; do
      if [[ "$type" == "AMR" ]]; then
        gset["$symbol"]=1
      fi
    done < <(tail -n +2 "$GENOME_FILE")
  fi
  if [[ -f "$PLASMID_FILE" ]]; then
    while IFS=$'\t' read -r _ contig start stop strand symbol _ _ type _; do
      if [[ "$type" == "AMR" ]]; then
        pset["$symbol"]=1
      fi
    done < <(tail -n +2 "$PLASMID_FILE")
  fi
  line="$sample"
  for sym in "${gene_symbols[@]}"; do
    if [[ -n "${gset[$sym]:-}" && -n "${pset[$sym]:-}" ]]; then
      val="GP"
    elif [[ -n "${gset[$sym]:-}" ]]; then
      val="G"
    elif [[ -n "${pset[$sym]:-}" ]]; then
      val="P"
    else
      val="NA"
    fi
    line+=$'\t'$val
  done
  echo -e "$line" >> "$OUTPUT"
done < "$SAMPLE_LIST"

echo "Finished writing combined AMR summary to $OUTPUT"
