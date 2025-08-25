#!/usr/bin/env bash

while read -r sampleid; 
do
    samplename=$(cut -f 1 <<< "$sampleid")
    bc=$(cut -f 2 <<< "$sampleid")
    loc=$(cut -f 3 <<< "$sampleid")
    echo "${samplename}"
    dorado aligner \
        --threads 16 --verbose \
        blatem.fasta \
        "${loc}/SQK-NBD114-24_${bc}.bam" > "bam/${samplename}_blatem_aligned.bam"

    samtools sort -@ 16 -o "bam/${samplename}_blatem_sorted.bam" "bam/${samplename}_blatem_aligned.bam"

    samtools index -@ 16 "bam/${samplename}_blatem_sorted.bam"

    samtools coverage -d 0 -o "cov/${samplename}_cov.txt" "bam/${samplename}_blatem_sorted.bam"

    samtools coverage -d 0 -m -D "bam/${samplename}_blatem_sorted.bam"

    bcftools mpileup \
        --config ont-sup-1.20 \
        --threads 16 -d 5000 --output-type v \
        --max-idepth 5000 \
        --annotate INFO/AD \
        -o "vcf/${samplename}_raw.vcf" \
        -f blatem.fasta \
        "bam/${samplename}_blatem_sorted.bam"
     
    bcftools call \
        --threads 16 \
        -A -m --ploidy 1 \
        -o "vcf/${samplename}_called.vcf" \
        "vcf/${samplename}_raw.vcf"

    bcftools filter \
        --output-type v --threads 16 \
        -i 'QUAL>=25' \
        "vcf/${samplename}_called.vcf" \
        | awk '$5!="."' > "vcf/${samplename}_final.vcf"

    echo

done < blatem_sample_id_bc.txt
