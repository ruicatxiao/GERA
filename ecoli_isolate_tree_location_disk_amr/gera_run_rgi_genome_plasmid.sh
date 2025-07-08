#!/usr/bin/env bash
# Running RGI 6.0.4 with latest CARD database via docker

# Pulling latest RGI with built in database
docker pull quay.io/biocontainers/rgi:6.0.4--pyh05cac1d_0

while read -r sampleid; 
do

    echo "${sampleid}"
    if test -f "genome/${sampleid}_genome.fasta"; then
        echo "${sampleid} genome fasta files exists."
        docker run -v $PWD:/data quay.io/biocontainers/rgi:6.0.4--pyh05cac1d_0 \
            rgi main \
                --input_sequence "/data/genome/${sampleid}_genome.fasta" \
                --output_file "/data/rgi/genome/rgi_genome_${sampleid}" \
                --clean \
                -n 8
    fi

    if test -f "plasmid/${sampleid}_plasmid.fasta"; then
        echo "${sampleid} plasmid fasta files exists."
        docker run -v $PWD:/data quay.io/biocontainers/rgi:6.0.4--pyh05cac1d_0 \
            rgi main \
                --input_sequence "/data/plasmid/${sampleid}_plasmid.fasta" \
                --output_file "/data/rgi/plasmid/rgi_plasmid_${sampleid}" \
                --clean \
                -n 8
    fi
    echo "${sampleid} is finished"
    echo
done < rig_card_sample.txt



docker run -v $PWD:/data quay.io/biocontainers/rgi:6.0.4--pyh05cac1d_0 \
    rgi heatmap \
    -i /data/rgi/genome \
    -o /data/rgi/genome/all_genome_rgi_card

docker run -v $PWD:/data quay.io/biocontainers/rgi:6.0.4--pyh05cac1d_0 \
    rgi heatmap \
    -i /data/rgi/plasmid \
    -o /data/rgi/plasmid/all_plasmid_rgi_card



