#!/bin/bash

# This is for running mob suite 

# mob suite version is 3.1.8

singularity pull mob_suite.sif https://depot.galaxyproject.org/singularity/mob_suite:3.1.8--pyhdfd78af_1

singularity exec \
 --bind "/data/lmattei/galapagos/analysis/mobtyper_isolate_plasmids:/data/input" \
 --bind "/data/lmattei/galapagos/analysis/mobtyper_isolate_plasmids:/data/output" \
 --cpus 16 \
 --containall \
 --workdir /tmp \
 --home /data/input \
 --memory 200G \
 /usr/local/bin/mob_suite.sif \
 mob_typer --multi \
      --infile /data/input/plasmids.fasta \
      --out_file /data/output/plasmids_mobtyper_results.txt