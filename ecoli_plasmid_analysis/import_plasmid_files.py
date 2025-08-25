import os
from Bio import SeqIO

input_directory = '/Users/arnav/Galapagos_plasmids_analysis/genomes_all/'
output_filename = '/Users/arnav/Galapagos_plasmids_analysis/Combined_galap_plasmids.fna'

with open(output_filename, "w") as outfile:
    for filename in os.listdir(input_directory):
       if filename.endswith(".fasta"):
           filepath = os.path.join(input_directory, filename)
           for record in SeqIO.parse(filepath, "fasta"):
               new_id = f"{filename[:-4]}_{record.id}"
               record.id = new_id
               SeqIO.write(record, outfile, "fasta")
