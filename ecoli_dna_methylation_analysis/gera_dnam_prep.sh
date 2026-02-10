#!/bin/bash

# Define input files and output names
GFF1="SW071620.gff"
GFF2="SW070910.gff"
BED_THREADS=20

# Step 1: Extract mRNA names from both GFFs
awk '$3=="mRNA"' "$GFF1" | grep "Name=" | cut -f 9 | awk -F ';' '{print $2}' | sed 's/Name=//g' | sort -V | uniq -c | awk '{print $2}' > "${GFF1%.gff}_mrna_name.txt}"
awk '$3=="mRNA"' "$GFF2" | grep "Name=" | cut -f 9 | awk -F ';' '{print $2}' | sed 's/Name=//g' | sort -V | uniq -c | awk '{print $2}' > "${GFF2%.gff}_mrna_name.txt}"

# Step 2: Find shared mRNA names across both genomes
cat "${GFF1%.gff}_mrna_name.txt" "${GFF2%.gff}_mrna_name.txt" \
    | sort -V | uniq -c | awk '$1==2' | awk '{print $2}' | sort -V > shared_gene_names.txt

# Step 3: Create BED files for gene bodies (GB) and promoters (50bp upstream)
for GFF in "$GFF1" "$GFF2"; do
    ID=$(basename "$GFF" .gff)
    echo "Processing $ID..."

    # Extract mRNA features with Name= and create GB BED
    grep -wFf shared_gene_names.txt "$GFF" | awk '$3=="mRNA"' | cut -f 1,4,5,7,9 \
        | awk '{split($5,a,";"); print $1,$2,$3,a[2],".",$4}' OFS='\t' \
        | sed 's/Name=//g' | awk '$4!="repB"' | sort -k1,1 -k2,2n > "${ID}_gb.bed"

    # Create promoter BED (50bp upstream)
    grep -wFf shared_gene_names.txt "$GFF" | awk '$3=="mRNA"' | cut -f 1,4,5,7,9 \
        | awk '{split($5,a,";"); print $1,$2,$3,a[2],".",$4}' OFS='\t' \
        | sed 's/Name=//g' | sort -k1,1 -k2,2n \
        | awk '{if ($6 == "+") print $1, ($2-50),$2,$4,$5,$6; else if($6 == "-") print $1,$3,($3+50),$4,$5,$6}' OFS='\t' \
        | awk '$4!="repB"' | sort -k1,1 -k2,2n > "${ID}_promoter.bed"

    # Optional: Check for invalid coordinates
    awk '$2<=0 || $3<=0' "${ID}_promoter.bed" > "${ID}_promoter_invalid.bed"
done

# Step 4: Convert BED files to bgzipped and indexed (for modkit)
for GFF in "$GFF1" "$GFF2"; do
    ID=$(basename "$GFF" .gff)
    echo "Compressing $ID..."

    # Check if the corresponding bed file exists
    if [ -f "${ID}_gb.bed" ]; then
        bgzip -k -@ "$BED_THREADS" "${ID}_gb.bed"
    fi

    if [ -f "${ID}_promoter.bed" ]; then
        bgzip -k -@ "$BED_THREADS" "${ID}_promoter.bed"
    fi

    # Index the bgzipped files
    tabix -p bed "${ID}_gb.bed.gz"
    tabix -p bed "${ID}_promoter.bed.gz"
done

# Step 5: Generate modkit stats for both genomes (GB and promoter regions)
for GFF in "$GFF1" "$GFF2"; do
    ID=$(basename "$GFF" .gff)

    # For gene body (gb)
    modkit stats --mod-codes "m" --min-coverage 20 --threads "$BED_THREADS" \
        --no-header --force --regions "${ID}_gb.bed" --out-table "${ID}_gb_5mC.tsv}" "$ID"_gb.bed.gz

    modkit stats --mod-codes "a" --min-coverage 20 --threads "$BED_THREADS" \
        --no-header --force --regions "${ID}_gb.bed" --out-table "${ID}_gb_6mA.tsv}" "$ID"_gb.bed.gz

    modkit stats --mod-codes "21839" --min-coverage 20 --threads "$BED_THREADS" \
        --no-header --force --regions "${ID}_gb.bed" --out-table "${ID}_gb_4mC.tsv}" "$ID"_gb.bed.gz

    # For promoter
    modkit stats --mod-codes "m" --min-coverage 20 --threads "$BED_THREADS" \
        --no-header --force --regions "${ID}_promoter.bed" --out-table "${ID}_promoter_5mC.tsv}" "$ID"_promoter.bed.gz

    modkit stats --mod-codes "a" --min-coverage 20 --threads "$BED_THREADS" \
        --no-header --force --regions "${ID}_promoter.bed" --out-table "${ID}_promoter_6mA.tsv}" "$ID"_promoter.bed.gz

    modkit stats --mod-codes "21839" --min-coverage 20 --threads "$BED_THREADS" \
        --no-header --force --regions "${ID}_promoter.bed" --out-table "${ID}_promoter_4mC.tsv}" "$ID"_promoter.bed.gz
done

# Cleanup: Optional – remove intermediate files (comment out if needed)
# rm *.txt *.bed *.bed.gz *.bed.gz.*

echo "Script completed successfully!"