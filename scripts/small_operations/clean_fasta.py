from Bio import SeqIO

genome = '../../data/genomic/pfs1_scaffolds_pb.fasta'
valid_contigs = list()

for record in SeqIO.parse(genome, "fasta"):
    length = len(record)
    if length > 1000:
        valid_contigs.append(record)

SeqIO.write(valid_contigs, '../../data/genomic/scaffolds_1000.fasta', 'fasta')
