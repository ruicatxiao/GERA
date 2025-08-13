#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import collections
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
import os

sensitive = 2000  # granularity of plots
output_dir = "heatmap_plots"
os.makedirs(output_dir, exist_ok=True)

col_names = ["qseqid", "qlen", "sseqid", "slen", "length", "pident", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "qcovs"]
Blast_results = pd.read_csv('Galap_Plasmid_Blast_Results.txt', sep="\t", low_memory=False, names=col_names, header=None)
AMR_ncbi = pd.read_csv('output_ncbi.txt', sep="\t", low_memory=False)
sequences = [">SW070910_reoriented.f_2", ">SW070912_reoriented.f_2", ">SW071620_reoriented.f_2", ">PCP06238_reoriented.f_2"]


index = 0

for seq in sequences:
    seq_name = seq.replace(">", "").replace("/", "_").replace(".", "_")
    
    AMR = AMR_ncbi.loc[AMR_ncbi["SEQUENCE"] == seq]
    Blast_results_one_genome = Blast_results[Blast_results["qseqid"] == seq]
    
    old_arange = np.arange(1, np.unique(Blast_results_one_genome["qlen"])[0] + 1)
    new_arange = np.arange(1, np.unique(Blast_results_one_genome["qlen"])[0] + 1, (np.unique(Blast_results_one_genome["qlen"])[0] + 1) / sensitive).astype(int)
    
    jaccard_array = np.zeros((len(new_arange), len(new_arange)))
    sets = collections.defaultdict(list)
    
    for i in tqdm(new_arange):
        test1 = Blast_results_one_genome["qstart"] < i + 0.1
        test2 = Blast_results_one_genome["qend"] > i - 0.1
        test3 = test1.astype(int) + test2.astype(int)
        sets[i].append(set(Blast_results_one_genome[test3 == 2]["sseqid"]))
    
    index_i = 0
    for i in tqdm(new_arange):
        index_j = 0
        for j in new_arange:
            if i > j - 1:
                if len(sets[i][0]) == 0 or len(sets[j][0]) == 0:
                    jaccard_array[index_i, index_j] = np.nan
                    jaccard_array[index_j, index_i] = np.nan
                else:
                    jaccard_calc = len(sets[i][0] & sets[j][0]) / len(sets[i][0] | sets[j][0])
                    jaccard_array[index_i, index_j] = jaccard_calc
                    jaccard_array[index_j, index_i] = jaccard_calc
            index_j += 1
        index_i += 1
    
    scaling_factor = np.unique(Blast_results_one_genome["qlen"])[0] / sensitive
    
    # original heatmap
    plt.close()
    plt.figure(index, figsize=(24, 24))  # Increased figure size
    fig, ax = plt.subplots(1, figsize=(24, 24))  # Increased figure size
    image = ax.imshow(jaccard_array, cmap="Reds")
    ax.set_xticks(np.linspace(0, sensitive, 10), np.linspace(0, np.unique(Blast_results_one_genome["qlen"])[0], 10).astype(int))
    ax.set_yticks(np.linspace(0, sensitive, 10), np.linspace(0, np.unique(Blast_results_one_genome["qlen"])[0], 10).astype(int))
    ax.tick_params(axis='x',rotation=90)
    fig.colorbar(image)
    
    for i in AMR.index:
        gene = AMR.loc[i]
        sta, sto = min(gene["START"], gene["END"]), max(gene["START"], gene["END"])
        rect = patches.Rectangle((sta / scaling_factor, 0), (sto - sta) / scaling_factor, sensitive, linewidth=1, edgecolor='navy', facecolor='none')
        ax.add_patch(rect)
        rect = patches.Rectangle((0, sta / scaling_factor), sensitive, (sto - sta) / scaling_factor, linewidth=1, edgecolor='navy', facecolor='none')
        ax.add_patch(rect)
        rect = patches.Rectangle((sta / scaling_factor, sta / scaling_factor), (sto - sta) / scaling_factor, (sto - sta) / scaling_factor, linewidth=1, edgecolor='navy', facecolor='navy')
        ax.add_patch(rect)
    
    square_svg_path = os.path.join(output_dir, f"square_heatmap_{seq_name}.svg")
    plt.savefig(square_svg_path, format='svg', bbox_inches='tight')
    print(f"Saved square heatmap to: {square_svg_path}")
    plt.show()
    
    # triangle Heatmap
    plt.figure(index + 100, figsize=(24, 24))
    fig, ax = plt.subplots(1, figsize=(24, 24))

    jaccard_array_triangle = jaccard_array.copy()
    mask = np.triu(np.ones_like(jaccard_array_triangle, dtype=bool), k=1)
    jaccard_array_triangle[mask] = np.nan

    image = ax.imshow(jaccard_array_triangle, cmap="Reds")
    ax.set_xticks(np.linspace(0, sensitive, 10), np.linspace(0, np.unique(Blast_results_one_genome["qlen"])[0], 10).astype(int))
    ax.set_yticks(np.linspace(0, sensitive, 10), np.linspace(0, np.unique(Blast_results_one_genome["qlen"])[0], 10).astype(int))
    ax.tick_params(axis='x',rotation=90)
    fig.colorbar(image)
    vertices = [(0, 0), (0, sensitive), (sensitive, sensitive)]
    clip_poly = patches.Polygon(vertices, closed=True, transform=ax.transData)

    for i in AMR.index:
        gene = AMR.loc[i]
        sta, sto = min(gene["START"], gene["END"]), max(gene["START"], gene["END"])
        
        rect3 = patches.Rectangle((sta / scaling_factor, sta / scaling_factor), (sto - sta) / scaling_factor, (sto - sta) / scaling_factor, linewidth=1, edgecolor='navy', facecolor='navy')
        rect3.set_clip_path(clip_poly)
        ax.add_patch(rect3)
        
        gene_name = gene["GENE"]
        mid_pos = (sta + sto) / 2 / scaling_factor
        
        ax.text(mid_pos, mid_pos - 75, gene_name, rotation=0, fontsize=12, color='black', 
                ha='center', va='top', 
                bbox=dict(boxstyle="round,pad=0.3", fc=(1, 1, 0, 0.5), ec="none"))

    triangle_svg_path = os.path.join(output_dir, f"triangle_heatmap_{seq_name}.svg")
    plt.savefig(triangle_svg_path, format='svg', bbox_inches='tight')
    print(f"Saved triangle heatmap to: {triangle_svg_path}")
    plt.show()
    
    index += 1