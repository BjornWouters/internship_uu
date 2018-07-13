# Gene coverage exploration script
# Developer: Björn Wouters
# Date: 26-6-2018
# Function: Describes the contigs that have abnormal coverage 

# load libraries
library(superheat)

# Thresholds
lowest_pfs1 = 0.01
lowest_pfs = 0.1
sliding_window = 5

# Set the working directory. 
setwd('Documents/school/internship/')

# Load the file 
cov_file <- read.table(file = 'data/gene_cov/Pfs_genecov_all.tsv', sep = '\t', header = TRUE)
cov_file <- cov_file[, order(names(cov_file))]
row.names(cov_file) <- cov_file$genename
cov_file <- cov_file[,-4]

# Testing:
# Only select specific columns 
cov_file <- cov_file[, c('Pfs1RZ', 'Pfs2RZ', 'Pfs3RZ_C', 'Pfs4PV_C', 'Pfs5RZ_C', 'Pfs6RZ', 'Pfs7PV', 
                         'Pfs8RZ_C', 'Pfs9RZ', 'Pfs10RZ', 'Pfs11PV', 'Pfs12RZ', 'Pfs13PV', 'Pfs14RZ', 
                         'Pfs15RZ', 'Pfs16PV_C')]

# Prepare the file 
# Remove the 0 values in the PFS1 column 
cov_file <- cov_file[!(cov_file$Pfs1RZ < lowest_pfs1),]
# Set all the under the threshold values to 0 and above the threshold value to 1
cov_file[cov_file < lowest_pfs] <- 0
cov_file[cov_file > lowest_pfs] <- 1

# Check for the reliability 

# Column variation 
plot(sapply(cov_file, var), xlab = 'Isolate', ylab = 'Variation')
text(sapply(cov_file, var), labels = colnames(cov_file), pos = 3)

# Column standard deviation
plot(sapply(cov_file, sd), xlab = 'Isolate', ylab = 'Standard deviation')
text(sapply(cov_file, sd), labels = colnames(cov_file), pos = 3)

# Column mean
plot(sapply(cov_file, mean), xlab = 'Isolate', ylab = 'Standard deviation')
text(sapply(cov_file, mean), labels = colnames(cov_file), pos = 3)

# Row variation 
median(sapply(t(cov_file), var))

###############
# Analysis 
###############

# Mutation matrix initialization 
resistoflay <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
califlay <- c(1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1)
clermont <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1)
campania <- c(1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0)
boeing <- c(1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0)
lazio <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0)
whale <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1)
pigeon <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0)
resistance_matrix <- as.data.frame(cbind(resistoflay, califlay, clermont, 
                                         campania, boeing, lazio, whale, pigeon))
gene_matrix <- matrix(1)

# Retrieve all the genes matching the right pattern in the specified window
for (species in colnames(resistance_matrix)) {
  mutation_vector <- resistance_matrix[species]
  column_names <- colnames(cov_file)

  if (is.null(colnames(gene_matrix))) {
    colnames(gene_matrix) <- species
  } else {
    gene_matrix <- cbind(gene_matrix, 0)
    colnames(gene_matrix)[ncol(gene_matrix)] <- species  
  }
  
  for (i in 1:as.integer(length(t(mutation_vector))-sliding_window+1)) {
    range <- i:as.integer(i+sliding_window-1)
    window <- mutation_vector[range,]
    if (length(subset(window, window == 0)) >= 2 & length(subset(window, window == 1)) >= 2) {
      gene_subset <- subset(cov_file, cov_file[column_names[range[1]]] == window[1] &
                              cov_file[column_names[range[2]]] == window[2] & 
                              cov_file[column_names[range[3]]] == window[3] & 
                              cov_file[column_names[range[4]]] == window[4] & 
                              cov_file[column_names[range[5]]] == window[5]) #&
                              #cov_file[column_names[range[6]]] == window[6])
      if (length(gene_subset) > 0) {
        for (name in rownames(gene_subset)) {
          if (name %in% rownames(gene_matrix)) {
            gene_matrix[,species][name] <- gene_matrix[,species][name] + 1
          } else {
            if (is.null(rownames(gene_matrix))) {
              rownames(gene_matrix) <- name
            } else {
              gene_matrix <- rbind(gene_matrix, 0)
              rownames(gene_matrix)[nrow(gene_matrix)] <- name  
              gene_matrix[,species][name] <- gene_matrix[,species][name] + 1
            }
          }
        }
      }
    }
  }
}
# Normalize to 0-100%
gene_matrix <- gene_matrix/11*100
gene_matrix[gene_matrix == 0] <- NA

superheat(gene_matrix, scale=TRUE, heat.col.scheme = "red", 
          heat.na.col = "white", left.label.size = 0.4, left.label.text.alignment = "left")




