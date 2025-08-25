#!/usr/bin/env bash

python gera_cov_mod.py \
    --fa_input ${genome_fa_input} \
    --output cov_mod \
    --sample_barcode sample_barcode_folder.txt \
    --threads 36