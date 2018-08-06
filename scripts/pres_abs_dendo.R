
# load libraries
library(pvclust)

# biocLite("ComplexHeatmap")

# Set the working directory. 
setwd('Documents/school/internship/repo/internship_uu/data/')
cov_file_pfs <- read.table('prepared_dataset.csv', sep = ',', row.names = 'mrna_id', header = 1)
cov_file_non_pfs <- read.table('prepared_test_dataset.csv', sep = ',', row.names = 'mrna_id', header = 1)
merged_file <- merge(cov_file_pfs, cov_file_non_pfs, by = "row.names", all = TRUE)
rownames(merged_file) <- merged_file[,'Row.names']
merged_file[,'Row.names'] <- NULL
merged_file[,'ES1513'] <- NULL

#full data range!!
# bootstraps the analysis to determine significance of observed cluster.
genecov.pv <- pvclust(merged_file, nboot=100, method.dist =
                        "euclidean", method.hclust = "ward.D2") 

#boston.pv <- pvclust(Boston, nboot=100) # bootstraps the analysis to determine significance of observed cluster.
plot(genecov.pv)

genecov_scaled = apply(cov_file, 2, scale) # columns are scaled to allow trends to be seen.
Heatmap(genecov_scaled, cluster_columns = genecov.pv$hclust,
        heatmap_legend_param = list(title = "Scaled Coverage")) ## this code combines heatmap and pvclust() generated dendrogram.
rownames(Genecov_final) <- c()
Heatmap(Genecov_final, cluster_columns = genecov.pv$hclust,
        heatmap_legend_param = list(title = "Coverage")) ## this code combines heatmap and pvclust() generated dendrogram.

#full binary range!!
binary.pv <- pvclust(binary, nboot=100) # bootstraps the analysis to determine significance of observed cluster.
#boston.pv <- pvclust(Boston, nboot=100) # bootstraps the analysis to determine significance of observed cluster.
plot(binary.pv)
#binary_scaled = apply(binary, 2, scale) # columns are scaled to allow trends to be seen.
Heatmap(binary_scaled, cluster_columns = binary.pv$hclust,
        heatmap_legend_param = list(title = "Scaled Binary Coverage")) ## this code combines heatmap and pvclust() generated dendrogram.
rownames(binary) <- c()
Heatmap(binary, cluster_columns = binary.pv$hclust,
        heatmap_legend_param = list(title = "Binary Coverage")) ## this code combines heatmap and pvclust() generated dendrogram.

#See wich contigs reside with missing genes!!!!!

