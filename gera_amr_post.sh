#!/bin/bash
set -euo pipefail

#############################################
# CONFIGURATION
#############################################
# Directory where the sample files are stored.
INPUTDIR="/data/ruicatxiao/gera_analysis/amr"  # <-- change this to your actual input directory

# File containing sample IDs (one per line)
SAMPLE_LIST="/data/ruicatxiao/gera_analysis/amr/amr_sampleid.txt"  # <-- update with your sample ID file path

# Output file names (defined at the top)
OUTPUT_GENOME="/data/ruicatxiao/gera_analysis/amr/genome_amr_summary.tsv"
OUTPUT_PLASMID="/data/ruicatxiao/gera_analysis/amr/plasmid_amr_summary.tsv"

# Array of AMR classes (in the desired order)
classes=("AMINOGLYCOSIDE" "BETA-LACTAM" "EFFLUX" "FOSFOMYCIN" "FOSMIDOMYCIN" "LINCOSAMIDE/MACROLIDE/STREPTOGRAMIN" "MACROLIDE" "PHENICOL" "QUINOLONE" "STREPTOTHRICIN" "SULFONAMIDE" "TETRACYCLINE" "TRIMETHOPRIM")

#############################################
# CREATE OUTPUT FILE HEADERS
#############################################
# Write header for genome summary file
{
  printf "SAMPLEID_GENOME"
  for cl in "${classes[@]}"; do
    printf "\t%s" "$cl"
  done
  printf "\n"
} > "$OUTPUT_GENOME"

# Write header for plasmid summary file
{
  printf "SAMPLEID_PLASMID"
  for cl in "${classes[@]}"; do
    printf "\t%s" "$cl"
  done
  printf "\n"
} > "$OUTPUT_PLASMID"

#############################################
# FUNCTION: Process a TSV file and count classes
#############################################
# This function uses awk to:
#   - Skip the header line (NR > 1)
#   - Only consider rows where the 9th column ("Type") equals "AMR"
#   - Check the 11th column ("Class") against the provided list of classes
#   - Increment a counter for each class match.
# The fields are assumed to be tab-delimited.
process_file() {
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

#############################################
# MAIN LOOP: Process Each Sample
#############################################
# Loop through each sample ID in the provided list
while IFS= read -r sample || [ -n "$sample" ]; do
    # Trim any leading/trailing whitespace from sample ID
    sample=$(echo "$sample" | xargs)

    #############################################
    # Process the Genome File
    #############################################
    GENOME_FILE="${INPUTDIR}/${sample}_genome_amrplus.tsv"
    if [[ -f "$GENOME_FILE" ]]; then
        # Process file and get class counts
        genome_result=$(process_file "$GENOME_FILE")
    else
        # If the file does not exist, output 13 zeros separated by tabs
        genome_result=$(printf "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0")
    fi
    # Append the sample ID and counts to the genome output file
    echo -e "${sample}\t${genome_result}" >> "$OUTPUT_GENOME"

    #############################################
    # Process the Plasmid File
    #############################################
    PLASMID_FILE="${INPUTDIR}/${sample}_plasmid_amrplus.tsv"
    if [[ -f "$PLASMID_FILE" ]]; then
        plasmid_result=$(process_file "$PLASMID_FILE")
    else
        plasmid_result=$(printf "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0")
    fi
    # Append the sample ID and counts to the plasmid output file
    echo -e "${sample}\t${plasmid_result}" >> "$OUTPUT_PLASMID"

done < "$SAMPLE_LIST"

echo "Processing complete. Results saved in $OUTPUT_GENOME and $OUTPUT_PLASMID."
