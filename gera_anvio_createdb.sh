#!/usr/bin/env bash

for genome in *.fasta; do
    base=$(basename "$genome" _genome.fasta)
    anvi-gen-contigs-database -f "$genome" -o "${base}.db" -T 1
    anvi-run-hmms -c "${base}.db"
done