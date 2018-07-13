# Gene coverage exploration script
# Developer: Björn Wouters
# Date: 10-7-2018
# Function: Explore the coverage of the contigs 

# Load libraries
library(seqinr)
library(stringi)
library(ggplot2) 

setwd('Documents/school/internship/')

contig_file <- read.fasta('data/pfs1_scaffolds_pb.fasta')

length <- unname(sapply(contig_file, length))
barplot(length)
abline(h=1000, col="red")
# 5029 contigs above 1k and 3606 lower than 1k basepears 
abline(v=5029, col="red")


# Check for genes on specific contig length 

gene_file <- read.fasta('data/Pfs1_Annotations.genomic.fasta')

contig_matrix <- matrix(sapply(contig_file, length))
# Add gene count column 
contig_matrix <- cbind(contig_matrix, 0)
colnames(contig_matrix) <- c('length', 'genes')

naming <- function(name) {
  (stri_extract(name, regex = 'Pb_[0-9]{1,4}'))
}

coverage <- function(name) {
  (as.numeric(gsub('cov_', '', stri_extract(name, regex = 'cov_[0-9]{1,3}\\.[0-9]*' ))))
}

contig_names <- sapply(names(contig_file), naming)
coverage_percentage <- sapply(names(contig_file), coverage)

rownames(contig_matrix) <- contig_names

for (gene in gene_file) {
  contig <- stri_extract(getAnnot(gene), regex = 'Pb_[0-9]{1,4}')
  contig_matrix[contig, ]['genes'] <- contig_matrix[contig, ]['genes'] + 1
}

contig_df <- as.data.frame(contig_matrix)
contig_df$names <- rownames(contig_df)
contig_df$dist <- contig_matrix[,'genes']/contig_matrix[,'length']*100

color.function <- colorRampPalette( c( "#CCCCCC" , "#C70039" ) )
color.ramp <- color.function( n = nrow( x = contig_df ) )
contig_df$color <-
  as.character(
    x = cut(
      x = rank( x = coverage_percentage )  # used to assign order in the event of ties
      , breaks = nrow( x = contig_df )  # same as the 'n' supplied in color.function()
      , labels = color.ramp  # label the groups with the color in color.ramp
    )
  )

barplot(contig_df$dist, names.arg = rownames(contig_matrix), col = contig_df$color, border=NA)
