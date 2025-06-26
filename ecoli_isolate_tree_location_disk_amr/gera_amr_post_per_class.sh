#!/bin/bash
set -euo pipefail




SAMPLE_LIST="gera_analysis/amr/amr_sampleid.txt"  # <-- update with your sample 
OUTPUT_GENOME="gera_analysis/amr/genome_amr_summary.tsv"
OUTPUT_PLASMID="gera_analysis/amr/plasmid_amr_summary.tsv"


classes=("AMINOGLYCOSIDE" "BETA-LACTAM" "EFFLUX" "FOSFOMYCIN" "FOSMIDOMYCIN" "LINCOSAMIDE/MACROLIDE/STREPTOGRAMIN" "MACROLIDE" "PHENICOL" "QUINOLONE" "STREPTOTHRICIN" "SULFONAMIDE" "TETRACYCLINE" "TRIMETHOPRIM")


  printf "SAMPLEID_GENOME"
  for cl in "${classes[@]}"; do
    printf "\t%s" "$cl"
  done
  printf "\n"
} > "$OUTPUT_GENOME"


{
  printf "SAMPLEID_PLASMID"
  for cl in "${classes[@]}"; do
    printf "\t%s" "$cl"
  done
  printf "\n"
} > "$OUTPUT_PLASMID"

    local file="$1"
    awk -F "\t" -v OFS="\t" '
    BEGIN {
        # Define the class order to be counted
        class[1] = "AMINOGLYCOSIDE"
        class[2] = "BETA-LACTAM"
        class[3] = "EFFLUX"
        class[4] = "FOSFOMYCIN"
        class[5] = "FOSMIDOMYCIN"
        class[6] = "LINCOSAMIDE/MACROLIDE/STREPTOGRAMIN"
        class[7] = "MACROLIDE"
        class[8] = "PHENICOL"
        class[9] = "QUINOLONE"
        class[10] = "STREPTOTHRICIN"
        class[11] = "SULFONAMIDE"
        class[12] = "TETRACYCLINE"
        class[13] = "TRIMETHOPRIM"
        # Initialize counters for each class
        for (i = 1; i <= 13; i++) {
            count[i] = 0
        }
    }
    NR > 1 {
        # Only process rows with Type equal to "AMR" (column 9)
        if ($9 == "AMR") {
            for (i = 1; i <= 13; i++) {
                if ($11 == class[i]) {
                    count[i]++
                }
            }
        }
    }
    END {
        # Output the counts for all classes in order, separated by tabs
        for (i = 1; i <= 13; i++) {
            printf "%s%s", count[i], (i < 13 ? OFS : "\n")
        }
    }' "$file"
}


while IFS= read -r sample || [ -n "$sample" ]; do

    sample=$(echo "$sample" | xargs)

    GENOME_FILE="${INPUTDIR}/${sample}_genome_amrplus.tsv"
    if [[ -f "$GENOME_FILE" ]]; then

        genome_result=$(process_file "$GENOME_FILE")
    else

        genome_result=$(printf "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0")
    fi

    echo -e "${sample}\t${genome_result}" >> "$OUTPUT_GENOME"


    PLASMID_FILE="${INPUTDIR}/${sample}_plasmid_amrplus.tsv"
    if [[ -f "$PLASMID_FILE" ]]; then
        plasmid_result=$(process_file "$PLASMID_FILE")
    else
        plasmid_result=$(printf "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0")
    fi

    echo -e "${sample}\t${plasmid_result}" >> "$OUTPUT_PLASMID"

done < "$SAMPLE_LIST"

echo "Processing complete. Results saved in $OUTPUT_GENOME and $OUTPUT_PLASMID."
