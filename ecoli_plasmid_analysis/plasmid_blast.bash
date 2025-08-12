makeblastdb -in /Users/arnav/40k_db.fna -parse_seqids -dbtype nucl

blastn -query /Users/arnav/Galapagos_plasmids_analysis/Combined_galap_plasmids.fna -db /Users/arnav/40k_db.fna -evalue 1e-200 -outfmt "6 qseqid qlen sseqid slen length pident qstart qend sstart send evalue bitscore qcovhsp qcovs" -out /Users/arnav/Galapagos/Galap_Plasmid_Blast_Results.txt
