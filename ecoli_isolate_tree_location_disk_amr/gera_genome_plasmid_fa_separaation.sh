#!/usr/bin/env bash

INDIR="gera_analysis/liftoff_fa_gff"
GENOMEDIR="gera_analysis/genome"
PLASMIDIR="gera_analysis/plasmid"


# extract genome fasta and their gff
for fa in $(find "$INDIR" -name '*\.fasta');
 do
    echo
    echo "$fa"
    faname=$(basename "$fa" .fasta)
    echo "start processing ${faname}"
    seqkit sort -l -r $fa | seqkit fx2tab -l -Q | awk -F '\t' -v a="$faname" '{print a,$1,$3,$2}' OFS='\t' | awk -F '\t' '{split($2,b," "); print $1,b[1],$3,$4}' OFS='\t' | awk '$3>4000000' | awk '{print $2,$4}' OFS='\t' | seqkit tab2fx > "${GENOMEDIR}/${faname}_genome.fasta"
    echo "genome fasta done"

    cat "${INDIR}/${faname}.gff" | awk '$1=="contig00001"' > "${GENOMEDIR}/${faname}_genome.gff"
    echo "genome gff done"
    echo
done

# extract plasmid fasta and their gff

for fa in $(find "$INDIR" -name '*\.fasta');
 do
    echo
    echo "$fa"
    faname=$(basename "$fa" .fasta)
    echo "start processing ${faname}"
    seqkit sort -l -r $fa | seqkit fx2tab -l -Q | awk -F '\t' -v a="$faname" '{print a,$1,$3,$2}' OFS='\t' | awk -F '\t' '{split($2,b," "); print $1,b[1],$3,$4}' OFS='\t' | awk '$3<=4000000' | awk '{print $2,$4}' OFS='\t' | seqkit tab2fx > "${PLASMIDIR}/${faname}_plasmid.fasta"
    echo "plasmid fasta done"

    cat "${INDIR}/${faname}.gff" | awk '$1!="contig00001"' > "${PLASMIDIR}/${faname}_plasmid.gff"
    echo "plasmid gff done"
    echo
done