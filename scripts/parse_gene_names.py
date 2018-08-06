import csv
import re

from Bio import SeqIO
import pandas as pd

annotation_file = '../data/annotation_information.tsv'
annotation_df = pd.read_csv(annotation_file, sep='\t', header=0, index_col='TranscriptId')

with open('../data/gene_identifiers.csv', 'w') as output:
    writer = csv.writer(output, delimiter=',')
    writer.writerow(['mrna_id', 'gene_name', 'contig_length', 'repeatregion'])
    for record in SeqIO.parse("../data/genomic/Pfs1_Annotations.genomic.fasta", "fasta"):
        contig_length = re.search('length_[0-9]*', record.description).group(0).strip('length_')

        repeat = False
        annotation_row = annotation_df.loc[record.id]

        if type(annotation_row) is not pd.Series:
            annotation = annotation_row.iloc[0].repeatregion
        else:
            annotation = annotation_row.repeatregion

        if type(annotation) is str:
            repeat = True

        writer.writerow([record.id, contig_length, repeat])
