#!/usr/bin/env bash
while read -r sampleid; 
do
    echo
    echo "${sampleid} step 1 start"

    seqkit fq2fa \
        "${sampleid}_calls_filtlong.fastq.gz" \
         -o "mags/${sampleid}/${sampleid}_reads.fasta"
    echo "${sampleid} step 1 done"

    echo
    echo "${sampleid} step 2 start"

    makeblastdb -dbtype nucl \
        -in "mags/${sampleid}/${sampleid}_reads.fasta" \
        -title "${sampleid}_db" \
        -out "mags/${sampleid}/${sampleid}_db"
    echo "${sampleid} step 2 done"
 
    echo
    echo "${sampleid} step 3 start"

    if test -f "mags/${sampleid}/${sampleid}_amr_prot.fasta"; then
      echo "${sampleid}_amr_prot file exists."
        tblastn \
          -query "mags/${sampleid}/${sampleid}_amr_prot.fasta" \
          -db "mags/${sampleid}/${sampleid}_db" \
          -outfmt 6 \
          -out "mags/${sampleid}/${sampleid}_amr_all.txt" \
          -evalue 1e-2 \
          -num_alignments 1000000000 -num_threads 4 \
          -qcov_hsp_perc 90
        cat "mags/${sampleid}/${sampleid}_amr_all.txt" \
            | cut -f 2 | sort -V | uniq -c \
            | awk '{print $2}' > "mags/${sampleid}/${sampleid}_amr_read_id.txt"
        seqkit grep -f "mags/${sampleid}/${sampleid}_amr_read_id.txt" \
            "mags/${sampleid}/${sampleid}_reads.fasta" > "mags/${sampleid}/${sampleid}_amr_reads.fasta"
    fi
    echo "${sampleid} step 3 done"
   
    echo
done < sample_id.txt
