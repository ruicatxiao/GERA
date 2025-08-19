#!/usr/bin/env python3

import os
import subprocess
from collections import defaultdict


AMR_CORD_FILE = "amr_cord.txt"
FASTA_FILES = [
    "PCP06238_4circos.fasta",
    "SW070910_4circos.fasta",
    "SW070912_4circos.fasta",
    "SW071620_4circos.fasta"
]
CIRCOS_LINKS_DIR = "circos_links"
ALL_LINKS_GREY = "all.links_grey.txt"
AMR_LINKS = "amr.links.txt"
AMR_ANNOTATION = "amr_annotation.txt"
GERA_LINKS = "gera_links.txt"


AMR_COLORS = {
    "aph(3'')-Ib": "255,0,0",      # Red
    "aph(6)-Id": "0,0,255",        # Blue
    "blaEC-19": "0,255,0",         # Green
    "blaTEM-1": "255,0,255",       # Magenta
    "dfrA14": "255,165,0",         # Orange
    "mph(A)": "128,0,128",         # Purple
    "sul2": "0,255,255",           # Cyan
    "tet(B)": "255,255,0"          # Yellow
}

# Genes to skip (appear on only one plasmid)
SKIP_GENES = {"aadA1", "dfrA1", "erm(B)", "sat2_fam"}


def process_amr():
    gene_positions = defaultdict(list)

    with open(AMR_CORD_FILE, 'r') as f:
        for line in f:
            plasmid, start, end, gene = line.strip().split('\t')
            if gene in SKIP_GENES:
                continue
            gene_positions[gene].append((plasmid, start, end))

    with open(AMR_LINKS, 'w') as link_file:
        for gene, positions in gene_positions.items():
            if len(positions) < 2:
                continue
            color = AMR_COLORS.get(gene, "128,128,128") + ",0.25"
            for i in range(len(positions)):
                for j in range(i + 1, len(positions)):
                    p1, s1, e1 = positions[i]
                    p2, s2, e2 = positions[j]
                    link_file.write(f"{p1}\t{s1}\t{e1}\t{p2}\t{s2}\t{e2}\tcolor={color}\n")

    with open(AMR_ANNOTATION, 'w') as annot_file:
        with open(AMR_CORD_FILE, 'r') as f:
            for line in f:
                plasmid, start, end, gene = line.strip().split('\t')
                if gene in SKIP_GENES:
                    continue
                annot_file.write(f"{plasmid}\t{start}\t{end}\tcircle_color=black thickness=20p\n")

def run_nucmer_alignments():
    os.makedirs(CIRCOS_LINKS_DIR, exist_ok=True)

    ids = [os.path.basename(f).replace("_4circos.fasta", "") for f in FASTA_FILES]

    with open(ALL_LINKS_GREY, 'w') as out_file:
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                id1 = ids[i]
                id2 = ids[j]
                prefix = os.path.join(CIRCOS_LINKS_DIR, f"{id1}_vs_{id2}")
                print(f"Processing {id1} vs {id2}...")

                # Run nucmer
                subprocess.run([
                    "nucmer",
                    "--prefix", prefix,
                    f"{id1}_4circos.fasta",
                    f"{id2}_4circos.fasta"
                ], check=True)

                delta_file = f"{prefix}.delta"
                filter_delta = f"{prefix}.filter.delta"
                coords_file = f"{prefix}.coords"

                # Filter delta
                with open(filter_delta, 'w') as fd:
                    subprocess.run(["delta-filter", "-i", "99", "-l", "1000", delta_file], stdout=fd, check=True)

                # Parse coords
                result = subprocess.run(["show-coords", "-rclTH", filter_delta], stdout=subprocess.PIPE, text=True, check=True)
                lines = result.stdout.strip().splitlines()

                for line in lines:
                    if line.startswith("=") or not line.strip():
                        continue
                    parts = line.split()
                    if len(parts) < 13:
                        continue
                    start1, end1 = int(parts[0]), int(parts[1])
                    start2, end2 = int(parts[2]), int(parts[3])
                    contig1, contig2 = parts[11], parts[12]

                    if start1 > end1:
                        start1, end1 = end1, start1
                    if start2 > end2:
                        start2, end2 = end2, start2

                    out_file.write(f"{contig1}\t{start1}\t{end1}\t{contig2}\t{start2}\t{end2}\tcolor=232,232,232,0.5\n")

                # Cleanup
                for f in [delta_file, filter_delta, coords_file]:
                    if os.path.exists(f):
                        os.remove(f)

def combine_files():
    with open(GERA_LINKS, 'w') as outfile:
        for fname in [ALL_LINKS_GREY, AMR_LINKS]:
            if os.path.exists(fname):
                with open(fname, 'r') as infile:
                    outfile.write(infile.read())

# === Main ===
if __name__ == "__main__":
    process_amr()
    run_nucmer_alignments()
    combine_files()
    print(f"Generated combined links file: {GERA_LINKS}")