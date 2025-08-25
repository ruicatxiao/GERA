#!/usr/bin/env bash

# making table to summerize all mod data from bedmethyl
for bedmethyl in $(find gera_analysis/cov_mod/mod -name '*\.bedmethyl');
 do
    echo
    name=$(basename "$bedmethyl" .bedmethyl)
    echo "$name"
        awk -v filename=${name}  '
        {
            contig = $1
            mod_type = $4
            nmod = $12
            ntotal = $10

            contigs[contig] = 1

            if (mod_type == "m") {
                m_mod[contig] += nmod
                m_total[contig] += ntotal
            } else if (mod_type == "a") {
                a_mod[contig] += nmod
                a_total[contig] += ntotal
            } else if (mod_type == "21839") {
                c_mod[contig] += nmod
                c_total[contig] += ntotal
            }
        }
        END {
            for (contig in contigs) {
                m_percent = (m_total[contig] ? (m_mod[contig] / m_total[contig]) * 100 : 0)
                a_percent = (a_total[contig] ? (a_mod[contig] / a_total[contig]) * 100 : 0)
                c_percent = (c_total[contig] ? (c_mod[contig] / c_total[contig]) * 100 : 0)
                printf "%s %.0f %.2f %.2f %.2f\n", filename, contig, m_percent, a_percent, c_percent
            }
        }' $bedmethyl | awk '{print $1"_"$2, $1, $2, $3,$4,$5}' OFS='\t' >> gera_analysis/cov_mod/report/gera_methyl_all.tsv
    echo
done


for fa in $(find /data/lmattei/galapagos/24_isolates/analysis/autocycler/genomes_all -name '*\.fasta');
  do
     echo
     faname=$(basename "$fa" .fasta)
     echo "$faname"
        cat "$fa" | seqkit fx2tab -l -g \
            | sort -k4,4rn \
            | awk -F '\t' '{print $1,$4,$5}' OFS='\t' \
            | awk -F '\t' '{if ($1 ~ "dnaA") print $1,$2,$3,"genome"; else if($1 ~ "repA") print $1,$2,$3,"plasmid"; else print $1,$2,$3,"leftover"}' OFS='\t' \
            | awk -F '\t' '{split($1,a," "); print a[1],$2,$3,$4}'  OFS='\t' \
            | awk -F '\t' -v var="$faname" '{print var, $1, $2, $3, $4, $5}' OFS='\t' \
            | awk -F '\t' '{split($1,b,"_"); print b[1],$2,$3,$4, $5}'  OFS='\t' \
            | awk -F '\t' '{print $1"_"$2,$1,$2,$3,$4,$5 }' OFS='\t' >> gera_analysis/cov_mod/report/gera_genome_plasmid_gc_length.tsv
done


# Making table for coverage info
for report in $(find gera_analysis/cov_mod/bam -name '*\.txt');
  do
     echo
     reportname=$(basename "$report" .txt)
     echo "$reportname"
     tail -n+2 "$report" | cut -f 1,7 \
     | awk -F '\t' -v var="$reportname" '{print var"_"$1, var, $1, $2}' OFS='\t' >> gera_analysis/cov_mod/report/gera_meandepth.tsv
done

# Sorting them
for tsv in $(find gera_analysis/cov_mod/report -name '*\.tsv');
  do
     echo
     tsvname=$(basename "$tsv" .tsv)
     sort -k1,1 -V "$tsv" > "gera_analysis/cov_mod/report/${tsvname}_sorted.tsv"
done

# merging of 3 separate reports
awk -F'\t' -v OFS='\t' '
    BEGIN {
        # Load genome data into an associative array
        file = "gera_genome_plasmid_gc_length_sorted.tsv";
        while ((getline < file) > 0) {
            key = $1 FS $2 FS $3;
            genome[key] = $4 FS $5 FS $6;  # Include the new classification column
        }
        close(file);
        
        # Load meandepth data into another associative array
        file = "gera_meandepth_sorted.tsv";
        while ((getline < file) > 0) {
            key = $1 FS $2 FS $3;
            meandepth[key] = $4;
        }
        close(file);
        
        # Print header
        print "ID\tSample\tContig\tm5C\tm6A\tm4C\tLength\tGC\tSeqtype\tMeanDepth";
    }
    {
        # For each line in the main file, build the key
        key = $1 FS $2 FS $3;
        
        # Retrieve genome and meandepth data, default to NA if missing
        gen_data = (key in genome) ? genome[key] : "NA\tNA\tNA";
        depth_data = (key in meandepth) ? meandepth[key] : "NA";
        
        # Print the merged line with all columns
        print $1, $2, $3, $4, $5, $6, gen_data, depth_data;
    }' gera_methyl_all_sorted.tsv > gera_genome_all_stats.tsv