#!/usr/bin/env bash

INPUTDIR="/data/ruicatxiao/gera_analysis/liftoff_fa_gff"
# remove exisiting file to avoid duplicate rows
rm genome_stats_4plotting.tsv

for fa in $(find "$INPUTDIR" -name '*\.fasta');
 do
    echo
    echo "$fa"
    faname=$(basename "$fa" .fasta)
    echo "start processing ${faname}"
    
    # seqkit stats -b -T $fa >> ../genome_stats.txt
    seqkit sort -l -r $fa | seqkit fx2tab -l -g | cut -f 1,4,5 | awk -F '\t' -v a="$faname" '{print a,$1,$2,$3}' OFS='\t' | awk -F '\t' '{split($2,b," "); print $1,b[1],$3,$4}' OFS='\t' | awk '{if ($3 > 3000000) print $1, $2,"genome",$3,$4; else print $1, $2,"plasmid",$3,$4}' OFS='\t' >> genome_stats_4plotting.tsv 

    echo "done"
    echo
done

R -e "library(ggplot2); library(dplyr); 
df <- read.table('genome_stats_4plotting.tsv', header=FALSE, sep='\t', 
       col.names=c('sample_name','contig_name','seq_type','seq_length','seq_gc')); 
df_genome <- df %>% filter(seq_type=='genome'); 
df_genome\$sample_name <- factor(df_genome\$sample_name, 
       levels = names(sort(tapply(df_genome\$seq_length, df_genome\$sample_name, sum), decreasing=TRUE))); 
p <- ggplot(df_genome, aes(x=sample_name, y=seq_length, fill=contig_name)) + 
     geom_bar(stat='identity', color='black') + 
     theme_minimal() + 
     theme(legend.position='none', axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(size = 8)) + 
     labs(x='Sample', y='Genome Length'); 
ggsave('genome_barchart.pdf', p)"

R -e "library(ggplot2); library(dplyr); 
df <- read.table('genome_stats_4plotting.tsv', header=FALSE, sep='\t', 
       col.names=c('sample_name','contig_name','seq_type','seq_length','seq_gc')); 
df_plasmid <- df %>% filter(seq_type=='plasmid'); 
totals <- tapply(df_plasmid\$seq_length, df_plasmid\$sample_name, sum); 
df_plasmid\$sample_name <- factor(df_plasmid\$sample_name, levels = names(sort(totals, decreasing=TRUE))); 
p <- ggplot(df_plasmid, aes(x=sample_name, y=seq_length, fill=contig_name)) + 
     geom_bar(stat='identity', position='stack', color='black') + 
     theme_minimal() + 
     theme(legend.position='none', axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(size = 8)) + 
     labs(x='Sample', y='Total Plasmid Length'); 
ggsave('plasmid_stacked_barchart.pdf', p)"





R -e "library(ggplot2); library(viridis); 
df <- read.table('genome_stats_4plotting.tsv', header=FALSE, sep='\t', 
       col.names=c('sample_name','contig_name','seq_type','seq_length','seq_gc')); 
p <- ggplot(df, aes(x=seq_gc, y=seq_length)) + 
  geom_density_2d() + 
  scale_x_continuous(limits=c(25,75)) + 
  scale_y_continuous(limits=c(0,6000000)) + 
  scale_fill_viridis(option='magma', discrete = TRUE) + 
  theme_minimal() + 
  labs(x='GC (%)', y='Sequence Length'); 
ggsave('2d_density_hexplot.pdf', p)"



R -e "library(ggplot2); 
df <- read.table('genome_stats_4plotting.tsv', header=FALSE, sep='\t', 
       col.names=c('sample_name','contig_name','seq_type','seq_length','seq_gc')); 
p <- ggplot(df, aes(x=seq_gc, y=seq_length, color=seq_type)) + 
  geom_point(size=3, alpha=0.8) + 
  scale_x_continuous(limits=c(25,75)) + 
  scale_y_continuous(limits=c(0,6000000)) + 
  scale_color_manual(values=c('genome'='#1f78b4','plasmid'='#E69F00')) + 
  theme_minimal() + labs(x='GC (%)', y='Sequence Length'); 
ggsave('xy_scatter_plot.pdf', p)"
