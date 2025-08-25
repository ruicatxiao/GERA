#!/usr/bin/env bash
set -euo pipefail
# anvio processing

conda activate anvio-7.1
for genome in *_genome.fasta; do
    base=$(basename "$genome" _genome.fasta)
    anvi-gen-contigs-database -f "$genome" -o "${base}.db" -T 8
    anvi-run-hmms -c "${base}.db"
done

# Compute completeness of genome
for db in *.db; do
    echo
    echo "$db"
    anvi-compute-completeness -c "$db" \
    --completeness-source Bacteria_71
    echo
done


anvi-gen-genomes-storage -e external-genome.txt \
                         -o ECOLI-GENOMES.db

anvi-pan-genome -g ECOLI-GENOMES.db \
                --project-name "GERA_Ecoli_Pan" \
                --num-threads 32 

anvi-compute-genome-similarity --external-genomes external-genome.txt \
                               --program pyANI \
                               --output-dir ANI \
                               --num-threads 32 \
                               --pan-db GERA_Ecoli_Pan/GERA_Ecoli_Pan-PAN.db
anvi-display-pan \
    -g ECOLI-GENOMES.db \
    -p GERA_Ecoli_Pan/GERA_Ecoli_Pan-PAN.db \
    --server-only

anvi-get-sequences-for-hmm-hits \
    --external-genomes external-genome.txt \
    --hmm-source Bacteria_71 \
    --list-available-gene-names

anvi-get-sequences-for-hmm-hits \
    --external-genomes external-genome.txt \
    -o concatenated-proteins.fa \
    --hmm-source Bacteria_71 \
    --gene-names ADK,AICARFT_IMPCHas,ATP-synt,ATP-synt_A,Adenylsucc_synt,Chorismate_synt,EF_TS,Exonuc_VII_L,GrpE,Ham1p_like,IPPT,OSCP,PGK,Pept_tRNA_hydro,RBFA,RNA_pol_L,RNA_pol_Rpb6,RRF,RecO_C,Ribonuclease_P,Ribosom_S12_S23,Ribosomal_L1,Ribosomal_L13,Ribosomal_L14,Ribosomal_L16,Ribosomal_L17,Ribosomal_L18p,Ribosomal_L19,Ribosomal_L2,Ribosomal_L20,Ribosomal_L21p,Ribosomal_L22,Ribosomal_L23,Ribosomal_L27,Ribosomal_L27A,Ribosomal_L28,Ribosomal_L29,Ribosomal_L3,Ribosomal_L32p,Ribosomal_L35p,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6,Ribosomal_L9_C,Ribosomal_S10,Ribosomal_S11,Ribosomal_S13,Ribosomal_S15,Ribosomal_S16,Ribosomal_S17,Ribosomal_S19,Ribosomal_S2,Ribosomal_S20p,Ribosomal_S3_C,Ribosomal_S6,Ribosomal_S7,Ribosomal_S8,Ribosomal_S9,RsfS,RuvX,SecE,SecG,SecY,SmpB,TsaE,UPF0054,YajC,eIF-1a,ribosomal_L24,tRNA-synt_1d,tRNA_m1G_MT \
    --return-best-hit \
    --get-aa-sequences \
    --concatenate

anvi-gen-phylogenomic-tree \
    -f concatenated-proteins.fa \
    -o phylogenomic-tree.txt

anvi-interactive \
    -t phylogenomic-tree.txt \
    -p phylogenomic-profile.db \
    --server-only \
    --manual \
    -A meta.txt