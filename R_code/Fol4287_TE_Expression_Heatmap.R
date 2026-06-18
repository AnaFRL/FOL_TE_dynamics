# Author: Ana Rodríguez López
# Date: 03/03/25
# TEcount results for RNAseq PDB wt 16h Microconidia, 60h Microconidia, 12h Mycelium

# Libraries
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(viridis)

# Read files 
counts <- read.table("pdbTC_RNAseq_TEcount_all.txt", header = T, sep = "\t")
te_length <- read.table("TEs_sequence_length_Foxy3.txt", header = T, sep = "\t")

# Create separate columns for TE, family and class
counts <- counts %>% 
  separate(gene.TE, into = c("TE_Gene", "family", "class"), sep = ":", remove = FALSE)

# Normalize read counts. ---- 
# Based on: https://www.mathworks.com/help/bioinfo/ug/identifying-differentially-expressed-genes-from-rna-seq-data.html
# Steps:
# a. Geometric mean: You establish a reference value for each gene/TE.
# b. Ratios: You calculate how the counts in each sample deviate in relation to this geometric mean.
# c. Size factors: You calculate an adjustment for each sample that compensates for differences in library sizes.
# d. Normalization: You divide the counts by the size factor for each sample, ensuring that differences between samples are not due to library size, but to the true relative expression of genes/TEs.

# 1. Calculate geometric means for each TE and gene across all samples
counts_geometric_mean <- counts %>%
  group_by(TE_Gene) %>%
  summarise(geometric_mean = exp(mean(log(count), na.rm = TRUE)))

# 2. Join geometric mean values to the counts dataframe and calculate ratio counts/geom.mean
counts_normalized <- counts %>%
  left_join(counts_geometric_mean, by = "TE_Gene") %>%  # left_join here will add geomtric mean values from counts_geometric_mean dataframe
  mutate(ratio = count / geometric_mean)  # calculate ratios (normalization by geometric mean)

# 3. Replace NaN and Inf values (produced by the 0/0 division) with 0 before continuing
counts_normalized <- counts_normalized %>%
  mutate(ratio = ifelse(!is.finite(ratio), 0, ratio)) # Replace NaN and Inf with 0 in ratio column from counts_normalized dataframe

# 4. Calculate size factor (median ratios for each sample)
size_factors <- counts_normalized %>%
  group_by(sample) %>%
  summarise(size_factor = median(ratio, na.rm = TRUE))

# 5. Join size factor values to the counts dataframe and calculate counts/size.factor
counts_final <- counts_normalized %>%
  left_join(size_factors, by = "sample") %>%
  mutate(normalized_count = count / size_factor)

# 6. Obtain a dataframe with TE data only. 
# Create a list with TE names
TE_list <- te_length %>% 
  select(TE) %>% 
  pull() # Extract column as a vector
# Filtrate dataframe for TEs
counts_TE <- counts_final %>%
  filter(TE_Gene %in% TE_list)

# 7. Normalized expression using TE length 
counts_TE <- counts_TE %>%
  left_join(te_length, by = c("TE_Gene" = "TE")) %>%  # Join counts_TE dataframe with te_length dataframe
  mutate(norm_by_length = (normalized_count / length) * 1000)  # Normalized by length. Multiply *1000 to convert values to RPK (reads per kilobase)

# 8. Add "Condition" column
counts_TE$Condition <- as.factor(gsub("R[0-9]+$", "", counts_TE$sample))

# HEATMAPS SORTED BY TE FAMILIES (CUT-OFF EXPRESSION, LOG TRANSFORMED DATA) ----

# Create a vector with all TE family ordering
TEs_order <- c("Copia", "Gypsy/Ty3", "Tad1","Foxy", "hAT", "Tc1/mariner", "Helitron", "Mule", "piggyBac", "Crypton", "Marsu")  

# Calculate the expression mean by group of TE and add the column to the dataframe
counts_TE <- counts_TE %>% 
  group_by(TE_Gene) %>% 
  mutate(mean_norm_by_length = mean(norm_by_length, na.rm = TRUE))

# Sort the data frame with TEs_order vector and sort by norm_by_length (descendent order)
counts_TE$family<- factor(counts_TE$family, levels = TEs_order) # Convert family column in a sort factor
counts_TE <- counts_TE %>% arrange(family, desc(mean_norm_by_length), TE_Gene) # Arrange the data frame by family and expression (descendent, from higher to lower), and finally by TE name (alphabetic).

# Apply a cut-off of 10 for the mean_norm_by_length column
counts_TE_filtered <- counts_TE %>% 
  filter(mean_norm_by_length >= 10)

# Create an expression matrix with the normalized values for the TEs 
mean_matrix <- counts_TE_filtered %>%
  select(TE_Gene, mean_norm_by_length) %>%  
  distinct() %>% 
  column_to_rownames(var = "TE_Gene") %>%  
  as.matrix()  # Convierte a matriz

# Use a logaritmic transformation
mean_matrix_log <- log2(mean_matrix + 1)  # +1 for avoiding log(0)
colnames(mean_matrix_log) <- "Mean" # Change column name

# Heatmap with values
max_val <- ceiling(max(mean_matrix_log, na.rm = TRUE)) # To determine the max value for the legend
TE_expression_cutoff <- pheatmap(mean_matrix_log, 
                                cluster_rows = FALSE,  
                                cluster_cols = FALSE,  
                                show_rownames = TRUE, 
                                show_colnames = TRUE, 
                                color = colorRampPalette(c("#410257", "#FFFFFF", "#FF8000"))(50),
                                legend_breaks = seq(0, 18, by = 3))
ggsave("TEexpression_meanConditions_cutoffexpression10_TEheatmap_family_sorted_logTransf_TCPDB.png", plot = TE_expression_cutoff, width = 2.7, height = 12, dpi = 320)


# Plot heatmap with RNAseq data of studied conditions, but removing the TEs under the cutoff 

# Create numeric matrix 
counts_matrix <- counts_TE_filtered %>%
  select(TE_Gene, sample, norm_by_length) %>%  
  pivot_wider(names_from = sample, values_from = norm_by_length) %>%
  column_to_rownames(var = "TE_Gene") %>% 
  as.matrix()

# Sort the columns
new_order <- c("W12MR1", "W12MR2", "W12MR3", "W16ER1", "W16ER2", "W16ER3", "W60ER1", "W60ER2", "W60ER3")  # Orden deseado
counts_matrix <- counts_matrix[, new_order]  

# Logaritmic transformed data
counts_matrix_sorted_log <- log2(counts_matrix + 1)  # +1 for avoiding log(0)

# Heatmap scaled by TE
scaled_byTE_plot_cutoff <- pheatmap(counts_matrix_sorted_log, 
                             cluster_rows = FALSE,  
                             cluster_cols = FALSE,  
                             scale = "row",  # Normalizing the matrix is done using the scale argument of the heatmap() function. It can be applied to row or to column. Here the row option is chosen, since we need to absorb the variation between row (TEs) But we cannot make the comparison between transposons. Here we can only see the differences in expression for a TE in the three different conditions studied. (I think it makes a Z-score by each row?)
                             show_rownames = TRUE, 
                             show_colnames = TRUE, 
                             color = colorRampPalette(c("#410257", "#FFFFFF", "#FF8000"))(50))
ggsave("scale_by_TE_heatmap_cutoffexpression10_family_sorted_logTransf_TCPDB.png", plot = scaled_byTE_plot, width = 5, height = 12, dpi = 320)
