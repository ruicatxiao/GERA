#!/bin/bash

# --- Configuration ---
# Pattern to match your files
PATTERN="SW*_*_*m*.tsv"
OUTPUT_FILE="combined_methylation_matrix.tsv"

# --- Script Start ---
set -e # Exit immediately if a command exits with a non-zero status

# Create a temporary directory for processing
TEMP_DIR=$(mktemp -d)
echo "Using temporary directory: $TEMP_DIR"

# 1. Identify files
FILES=($PATTERN)
if [ ${#FILES[@]} -eq 0 ]; then
    echo "Error: No files found matching pattern '$PATTERN'"
    exit 1
fi
echo "Found ${#FILES[@]} files to process."

# 2. Process each file: Extract Gene (col 4) and Percent (col 8), then Sort by Gene
# We save the processed files into the temp directory.
SORTED_FILES=()
HEADER="Gene"

for f in "${FILES[@]}"; do
    # Construct clean column name from filename
    COL_NAME="${f%.tsv}"
    HEADER="$HEADER"$'\t'"$COL_NAME"
    
    # Define temp output path
    OUT_PATH="$TEMP_DIR/${COL_NAME}.sorted"
    
    # Extract columns 4 and 8, then sort by column 1 (Gene Name)
    # LC_ALL=C ensures standard ASCII sorting
    awk -F'\t' '{print $4 "\t" $8}' "$f" | LC_ALL=C sort -k1,1 > "$OUT_PATH"
    
    SORTED_FILES+=("$OUT_PATH")
done

# 3. Validation: Check if Gene Lists are identical across all files
echo "Checking gene name consistency across files..."

# Get the gene list from the first processed file (column 1 only)
BASE_GENE_LIST="${TEMP_DIR}/baseline_genes.txt"
cut -f1 "${SORTED_FILES[0]}" > "$BASE_GENE_LIST"

# Compare base list with every other file
ALL_OK=true
for i in "${!SORTED_FILES[@]}"; do
    CURRENT_FILE="${SORTED_FILES[$i]}"
    CURRENT_GENES="${TEMP_DIR}/current_genes.txt"
    
    # Extract genes from current file
    cut -f1 "$CURRENT_FILE" > "$CURRENT_GENES"
    
    # Compare with baseline
    if ! diff -q "$BASE_GENE_LIST" "$CURRENT_GENES" >/dev/null 2>&1; then
        echo "ERROR: Gene list mismatch detected!"
        echo "Baseline file: ${FILES[0]}"
        echo "Mismatching file: ${FILES[$i]}"
        echo "Check for typos, duplicates, or missing entries in these files."
        ALL_OK=false
        break
    fi
done

if [ "$ALL_OK" = false ]; then
    echo "Validation failed. Cleaning up and exiting."
    rm -rf "$TEMP_DIR"
    exit 1
fi

echo "Validation successful. All gene lists are identical."

# 4. Merge the files
# We paste the Gene column from the first file, and the Methylation column (col 2) from ALL files.
echo "Merging files..."

# Construct the paste command
# Start with gene column from the first file
CMD="paste <(cut -f1 ${SORTED_FILES[0]})"

# Add methylation column (col 2) from every file
for f in "${SORTED_FILES[@]}"; do
    CMD="$CMD <(cut -f2 $f)"
done

# Execute paste command and add header
eval "$CMD" > "$TEMP_DIR/merged_data_no_header.txt"

{
    echo -e "$HEADER"
    cat "$TEMP_DIR/merged_data_no_header.txt"
} > "$OUTPUT_FILE"

# 5. Finalize
rm -rf "$TEMP_DIR"
echo "--------------------------------------------------"
echo "Done! Output saved to: $OUTPUT_FILE"
echo "Total genes (rows): $(($(wc -l < "$OUTPUT_FILE") - 1))"
echo "--------------------------------------------------"
