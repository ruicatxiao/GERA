# GERA
Scripts and analysis used for GERA project
![Image](https://github.com/user-attachments/assets/eb9b428c-5938-486d-900a-fb9b93c016a4)

## Abstract
Since the publication of Charles Darwins On the Origin of Species in 1859, the Galapagos archipelago has become emblematic of natural selection and evolution. While the lens of evolution in the Galapagos has traditionally focused on iconic megafauna, including finches, marine iguanas, and giant tortoises, the marine environment is also home to diverse microbial ecosystems that are constantly evolving under selective pressure from environmental factors such as human activity. We focused on the second most populated island within the archipelago - San Cristobal - an island that has experienced rapid urbanization and intense international tourism pressure. Using a lab-free approach, we spatiotemporally mapped wastewater contamination around San Cristobal. On-site metagenomic sequencing revealed a stark shift in genera and a higher count of antimicrobial resistance genes at wastewater-associated sites. Over 40% of lactose-fermenting Enterobacteriaceae isolates collected from sewage and wastewater outfall exhibited multidrug resistance (MDR). Long-read sequencing and de novo assembly of bacterial genomes and plasmids from MDR Escherichia coli revealed frequent and rapid reassortment of antimicrobial resistance genes on plasmids, generating unique antimicrobial resistance profiles. This study not only provides a framework for conducting antimicrobial resistance research in low-resource settings but also underscores the impact of wastewater contamination on the environmental AMR landscape and highlights potential threats to human and animal health.



## Per Folder Scripts and Files Overview
Detailed per folder script explanation.

### Ecoli_isolate_basecall_assembly.
This folder contains scripts for basecalling and assembly of long-read sequencing data from E. coli isolates.

    gera_dorado_demux_call.sh
    Runs Dorado basecaller on POD5 files, demultiplexes barcoded reads, and converts BAMs to FASTQ format.

    gera_dnaapler.sh
    Uses DNAapler to reorient and annotate contigs from long-read assemblies.

    gera_autocycler.sh
    Runs multiple assemblers on subsampled reads using Autocycler assembler for comparative genome assembly.

### Ecoli_isolate_tree_location_disk_amr.
Scripts used for analyzing E. coli isolate data including phylogenetic context, location, antimicrobial resistance (AMR), and genomic features.

    gera_run_anvio.sh
    Generates contigs databases and computes completeness for genomes using Anvi'o.

    gera_run_rgi_genome_plasmid.sh
    Uses RGI to screen genomes and plasmids for resistance genes and generates heatmaps.

    gera_amr_run.sh
    Runs AMRFinder+ on genome and plasmid FASTA files to identify AMR genes.

    gera_amr_post_per_gene.sh
    Post-processes AMRFinder+ output and merges gene-level results.

    gera_amr_post_per_class.sh
    Post-processes AMRFinder+ output and summarizes results by AMR gene class.

    gera_amr_plotting.R
    Generates AMR heatmaps using R and visualization libraries.

    gera_amr_gff_prep.sh
    Prepares GFF files for genome annotation for AMR analysis.

    gera_amr_extract_protein.sh
    Extracts protein sequences from GFF annotation.

    gera_cov_mod_post.sh
    Processes coverage and methylation data, generating tables and combined summaries.

    gera_genome_plasmid_fa_separaation.sh
    Separates genome and plasmid FASTA files from the analysis.

    gera_genome_stats_plotting.sh
    Extracts and plots genomic and plasmid statistics, including GC content, length, and sequence type.

    gera_run_cov_mod.sh
    Wrapper script to run ONT DNA methylation and coverage analysis using gera_cov_mod.py.

    gera_cov_mod.py
    Python script for running Dorado aligner, samtools, and modkit to extract methylation and coverage data.

    gera_liftoff.sh
    Uses Liftoff to lift over annotations from a reference genome to a query genome.

    gera_gff_cleaning.sh
    Cleans GFF files using AGAT for compatibility with downstream tools.

    gera_anvio_createdb.sh
    Generates contigs databases for Anvi'o from FASTA files.

### Ecoli_plasmid_analysis
Contains scripts and files used for analyzing plasmids, including comparisons, annotations, and visualization using Circos.

    plasmid_blast.bash
    Runs BLASTN to compare plasmids against a database and outputs alignment results.

    import_plasmid_files.py
    Combines multiple FASTA files from different plasmid isolates into a single FASTA file, renaming records.

    gera_circos_links_highlights_prep.py
    Prepares Circos input files including links, highlights, and annotations for plasmid visualization.

    gera_circos.conf
    Configuration file for Circos plotting, defines karyotype, ideograms, and visual styles.

    gera_links.txt
    Contains Circos link data for plotting sequence similarities between plasmids.

    plasmid_oriv.txt
    Defines oriv regions in plasmids for Circos visualization.

    amr_annotation.txt
    Annotations for AMR genes in plasmids for Circos highlights.

    gera_run_mob_suite.sh
    Runs the MOB suite to classify plasmids by mobility type using Singularity.

    galapagos_plasmid_analysis_GERA.py
    Analyzes plasmid data in a heatmap format using Python, including visualizing AMR regions and SNP clustering.

### Metagenomic
Scripts for processing metagenomic data, especially long-read assemblies and AMR gene identification. 

    run_metaflye.sh
    Runs MetaFlye assembler on long-read data to generate metagenomic assemblies.

    run_amr_reads_prep.sh
    Prepares AMR reads by converting FASTQ to FASTA, building BLAST databases, and extracting reads with AMR hits.


### SNP_analysis (deprecated)
Deprecated scripts for SNP calling and plotting, particularly from AMR-focused reads.

    gera_amr_snp_vcf.sh
    Performs SNP calling from AMR-focused reads using Dorado, samtools, and bcftools.

    snp_plotting.R
    Visualizes SNP distribution, frequency, and clustering using ggplot2 and pheatmap.


## Manuscript Authors
```
Arnav Lal, Jade C. Riopelle, Katherine Villarin, 
Maya Mathur, Lia Enriquez, Naomi Phemister-Jimenez, 
Katherine Gilbert, Rui Xiao, Demy Castillo, 
Stephen Cole, Maddie Tilyou, Kelly Kennedy, 
Ernesto Vaca, Fausto Rodriguez, Whitman Cox, 
Wilson Castillo, Michael Weisberg, 
Lisa M. Mattei, Daniel P. Beiting
```

## Code Contributors
- Arnav Lal (https://github.com/Arnavlal)
- Lisa Mattei (https://github.com/lisa-mattei)
- Daniel Beiting (https://github.com/dpbisme)
- Rui Xiao (https://github.com/ruicatxiao)

## Data availability
- The ONT metagenomic reads are available through PRJNA132068
- The ONT isolate genomic reads and genome assemblies are available through PRJNA1320685




## Citations
```
"Human wastewater contamination drives the emergence of novel multidrug resistant bacteria in the Galápagos marine ecosystem" Arnav Lal, Jade C. Riopelle, Katherine Villarin, Maya Mathur, Lia Enriquez, Rui Xiao, Naomi Phemister-Jimenez, Katherine Gilbert, Stephen D. Cole, Maddie Tilyou, Kelly P. Kennedy, Ernesto Vaca, Wilson Castillo, Michael Weisberg, Lisa M. Mattei, Daniel P. Beiting bioRxiv 2025.09.12.675863; doi: https://doi.org/10.1101/2025.09.12.675863
```
