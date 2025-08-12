# Plasmid Analysis

This is incomplete, will fill in more details later

## Analysis 1: Comparison of Galapagos Plasmids to Global database

- Edit plasmid sequences --> w/ import_plasmid_files.py
- Blast against Global database --> w/ plasmid_blast.bash 
- 40,757 plasmid database, source: https://www.biorxiv.org/content/10.1101/2024.09.03.610986v1
- Generate heatmaps --> w/ galapagos_plasmid_analysis_GERA.py, using Blast output and AMR analysis of plasmids (AMRfinderplus & RGI, DOUBLE CHECK).

## Analysis 2: Clinker synteny analysis of plasmids

- Run Bakta on plasmids (DOUBLE CHECK version)
- Keep .gbff files
- Keep just AMR genes. 
- Run clinker (v0.0.31) on whole sequence, and just AMR genes, e.g. ```  clinker /.../*.gbff -p visualization.html ```
