# CNV calling for chromosomes

# Author: MVAP
# Version: 2023-03-08
# Adapted by ARL (2025-01-20)

# Directory tree

#create_subdirectories()

# Packages ----
# Remember to insert %>% is ctrl + Shift +M

ibrary(tidyverse)
library(readr)
library(dplyr)
library(ggplot2)

# Read files ----

# Location of bedgraph. There should be one per sample
# if you have bam but not bedgraph run at the command line 
# for i in `ls $FSCRATCH/*.bam`; do echo $i; bedtools genomecov -d -ibam $i > ${i}.bedgraph; done

setwd("/data/ana")

# get all files ending in "bedgraph"
bed.list <- list.files(pattern = "bedgraph")

# split name accordingly. the names will be use later in the graph
# Example:
# B138_...bedgraph / YC1_...bedgraph
fileList <- tibble(
  file = bed.list, 
  sample = gsub("_.*", "", bed.list)) # Extract sample name

# Initialize empty dataframe to store coverage data
bed.table <- data.frame()

for (i in 1:nrow(fileList)){
  print(i) # progress indicator (useful for large files)
  hname <- fileList$sample[i] # short name sample for figure
  
  # Read bedgraph file
  bed.file <- read.table(fileList$file[i], header = FALSE, sep = "\t", 
                         stringsAsFactors = FALSE,quote = "",
                         col.names = c("Scaffold", "Bp", paste(hname, "cov", sep="_")))
  # Build combined table iteratively
  if (dim(bed.table)[1] == 0){
    bed.table <- rbind(bed.table, bed.file) 
  }else{
    bed.table <- bed.table %>% left_join(bed.file, by = c("Scaffold"="Scaffold", "Bp"="Bp"))
  }
} # The last part of the code (if (dim)...) checks whether bed.table is empty. If bed.table has no rows (first iteration or uninitialized), it initializes bed.table with the content of bed.file using rbind(). Otherwise, it merges (left_join) the new bed.file with the existing bed.table by matching the columns "Scaffold" and "Bp", adding coverage data for the current sample to the corresponding positions. This structure allows progressively building a combined table of coverage across all samples.

# Check for NA values before proceeding
# any(is.na(bed.table))
# which(is.na(bed.table), arr.ind = TRUE)

# Functions ----

# Divide genome into 10 kb windows
rep10<-function(d){
  d$tenkb =  rep(1:ceiling(nrow(d)/10000), each = 10000)[1:nrow(d)]
  d
}

# Define pseudochromosome mapping
chromosome_map <- tibble(
  Scaffold= c("scaffold_1", "scaffold_4", "scaffold_16", "scaffold_15", "scaffold_17", "scaffold_26", "scaffold_7", "scaffold_41", "scaffold_19", "scaffold_34", "scaffold_31", "scaffold_20", "scaffold_27", "scaffold_28", "scaffold_18", "scaffold_24", "scaffold_35", "scaffold_29", "scaffold_2", "scaffold_3", "scaffold_5", "scaffold_6", "scaffold_9", "scaffold_8", "scaffold_10", "scaffold_11", "scaffold_12", "scaffold_13", "scaffold_32", "scaffold_14", "scaffold_22", "scaffold_38", "scaffold_23", "scaffold_37", "scaffold_43", "scaffold_39", "scaffold_25", "scaffold_42", "scaffold_33", "scaffold_40", "scaffold_44", "scaffold_21", "scaffold_30", "scaffold_36"),
  Chromosome = c("Chr1","Chr2","Chr2","Chr3","Chr3","Chr3","Chr3","Chr3","Chr3","Chr3","Chr3","Chr3","Chr3","Chr3","Chr3","Chr3","Chr3","Chr3","Chr4","Chr5","Chr7","Chr8","Chr9","Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Chr15","Chr15","Chr15","Chr15","Chr15","Chr15","Chr15","Chr15","Chr15","Chr15","Chr15","Chr15","Unpositioned","Unpositioned","Unpositioned"))


# CNVs ANALYSIS - MEDIAN APPROACH (no-WT normalized data)----
  
#Create 10 kb windows
bed.file_chromosomes <- bed.table %>% group_by(Scaffold) %>% group_modify(~rep10(.)) %>% ungroup()

# Calculate median per window
bed.median_chromosomes <- bed.file_chromosomes %>% group_by(Scaffold, tenkb) %>%   summarise(across(ends_with("cov"), ~median(.x))) %>% ungroup()

# Calculate global median per sample
bed.median_chromosomes_global <- bed.table %>%
  summarise(across(ends_with("cov"), ~median(.x, na.rm = TRUE)))

# Normalize by global median
data_columns <- colnames(bed.median_chromosomes)[-c(1, 2)] # Exclude 'Scaffold' and 'tenkb'
bed.medianNorm_chromosomes <- bed.median_chromosomes %>%
  mutate(across(all_of(data_columns), 
                ~ . / (bed.median_chromosomes_global %>% pull(cur_column()))))

  ## Plotting data non-WT normalized ---- 

# Keep scaffold order within the pseudochromosome
chromosome_map <- chromosome_map %>%
  mutate(Scaffold = factor(Scaffold, levels = unique(Scaffold)))

# Create window coordinates
bed.size <- bed.file_chromosomes %>% 
  group_by(Scaffold, tenkb) %>% 
  summarise(end = max(Bp), .groups = "drop") %>% 
  mutate(start =(tenkb*10000)-9999)

# Compute scaffold offsets for pseudochromosome layout
scaffold_offsets <- bed.size %>%
  group_by(Scaffold) %>%
  summarise(scaffold_len = max(end), .groups = "drop") %>%
  left_join(chromosome_map, by = "Scaffold") %>%
  # Sort by chromosome and by Scaffold
  arrange(Chromosome, Scaffold) %>%
  group_by(Chromosome) %>%
  # Calculate the offset
  mutate(offset = lag(cumsum(scaffold_len), default = 0)) %>%
  ungroup()

# Build plotting dataframe
bed.plot <- bed.medianNorm_chromosomes %>%
  left_join(bed.size, by = c("Scaffold", "tenkb")) %>%
  left_join(chromosome_map, by = "Scaffold") %>%
  left_join(dplyr::select(scaffold_offsets, Scaffold, offset), by = "Scaffold") %>%
  mutate(
    chr_start = start + offset,
    chr_end = end + offset
  )

# Long format
bed.t <- bed.plot %>%
  gather("sample", "cov", -c(Chromosome, Scaffold, tenkb, start, end, offset, chr_start, chr_end))

# Factor ordering
bed.t$Chromosome <- factor(
  bed.t$Chromosome,
  levels = c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr7","Chr8","Chr9",
             "Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Unpositioned"))

bed.t$sample<- factor(
  bed.t$sample, 
  levels= c("SRR7690004_cov", "YC1_cov", "YC2_cov", "YC3_cov","YC4_cov","YC5_cov","YP1_cov","YP2_cov", "YP3_cov", "YP4_cov", "YP5_cov"))

# Plot
ggplot(data=bed.t, aes(x=chr_start/1000000, y=cov, group=sample, color=sample)) + 
  geom_col() +
  facet_grid(sample ~ Chromosome, scales="free_x", space="free_x", labeller=label_value) +  
  xlab("Position in the genome (Mb)") + 
  ylab("Coverage density") + 
  theme_bw() +
  ggtitle(label = "Fol4287 YP (mycelium passaged) abd YC (conidia-passaged) samples Genome coverage to 2025 Assembly (Ayhan et al.)",
          subtitle = "Median every 10kb per Scaffold, each sample normalized by its median") +
  theme(axis.text.x = element_text(size=6),
        strip.text.x = element_text(size = 8),
        axis.text.y = element_text(size=6)) +
  scale_x_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 1)) +
  scale_color_manual(values = c(
    "SRR7690004_cov" = "black",  
    "YC1_cov" = "#90BFF9", "YC2_cov" = "#90BFF9", "YC3_cov" = "#90BFF9", "YC4_cov" = "#90BFF9", "YC5_cov" = "#90BFF9",  
    "YP1_cov" = "#DE8BF9", "YP2_cov" = "#DE8BF9", "YP3_cov" = "#DE8BF9", "YP4_cov" = "#DE8BF9", "YP5_cov" = "#DE8BF9"   
  ))

ggsave("CNVS_median_MedianNorm_pseudochromosomes.png", width = 35, height = 30, units = "in", dpi = 600)


# CNVs ANALYSIS - MEDIAN APPROACH AND NORMALIZATION BY WT ----

# Create 10 kb windows
bed.file_chromosomes <- bed.table %>% group_by(Scaffold) %>% group_modify(~rep10(.)) %>% ungroup()

# Calculate median per window
bed.median_chromosomes <- bed.file_chromosomes %>% group_by(Scaffold, tenkb) %>%   summarise(across(ends_with("cov"), ~median(.x))) %>% ungroup()

# Calculate global median for each sample
bed.median_chromosomes_global <- bed.table %>%
  summarise(across(ends_with("cov"), ~median(.x, na.rm = TRUE)))

# Normalize by global median
data_columns <- colnames(bed.median_chromosomes)[-c(1, 2)] # Excluir 'Scaffold' y 'tenkb'
bed.medianNorm_chromosomes <- bed.median_chromosomes %>%
  mutate(across(all_of(data_columns), 
                ~ . / (bed.median_chromosomes_global %>% pull(cur_column()))))

# Normalize by WT
wt_column <- "SRR7690004_cov"
bed.medianNorm_chromosomes_WTNorm <- bed.medianNorm_chromosomes %>% 
  mutate(across(all_of(data_columns), ~ . / .data[[wt_column]]))
# Replace NaN and Inf values (produced by the 0/0 division) with 0 before continuing
bed.medianNorm_chromosomes_WTNorm <- bed.medianNorm_chromosomes_WTNorm %>%
  mutate(across(everything(), ~ ifelse(is.nan(.) | is.infinite(.), 0, .))) # Replace NaN and Inf with 0

# Check that median in WT column = 1
bed.medianNorm_chromosomes_WTNorm %>%
  summarise(across(all_of(wt_column), median))

#Plotting coverages normalized to WT (everything normalized to WT, so WT is 1)

# Keep scaffold order within the pseudochromosome
chromosome_map <- chromosome_map %>%
  mutate(Scaffold = factor(Scaffold, levels = unique(Scaffold)))

# Create window coordinates
bed.size <- bed.file_chromosomes %>% 
  group_by(Scaffold, tenkb) %>% 
  summarise(end = max(Bp), .groups = "drop") %>% 
  mutate(start =(tenkb*10000)-9999)

# Compute scaffold offsets for pseudochromosome layout
scaffold_offsets <- bed.size %>%
  group_by(Scaffold) %>%
  summarise(scaffold_len = max(end), .groups = "drop") %>%
  left_join(chromosome_map, by = "Scaffold") %>%
  # Sort by chromosome and by Scaffold
  arrange(Chromosome, Scaffold) %>%
  group_by(Chromosome) %>%
  # Calculate the offset
  mutate(offset = lag(cumsum(scaffold_len), default = 0)) %>%
  ungroup()

# Build plotting dataframe
bed.plot <- bed.medianNorm_chromosomes_WTNorm %>%
  left_join(bed.size, by = c("Scaffold", "tenkb")) %>%
  left_join(chromosome_map, by = "Scaffold") %>%
  left_join(dplyr::select(scaffold_offsets, Scaffold, offset), by = "Scaffold") %>%
  mutate(
    chr_start = start + offset,
    chr_end = end + offset
  )

# Long format
bed.t <- bed.plot %>%
  gather("sample", "cov", -c(Chromosome, Scaffold, tenkb, start, end, offset, chr_start, chr_end))

# Factor ordering
bed.t$Chromosome <- factor(
  bed.t$Chromosome,
  levels = c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr7","Chr8","Chr9",
             "Chr10","Chr11","Chr12","Chr13","Chr14","Chr15","Unpositioned"))

bed.t$sample<- factor(
  bed.t$sample, 
  levels= c("SRR7690004_cov", "YC1_cov", "YC2_cov", "YC3_cov","YC4_cov","YC5_cov","YP1_cov","YP2_cov", "YP3_cov", "YP4_cov", "YP5_cov"))

# Plot
ggplot(data=bed.t, aes(x=chr_start/1000000, y=cov, group=sample, color=sample)) + 
  geom_col() +
  facet_grid(sample ~ Chromosome, scales="free_x", space="free_x", labeller=label_value) +  
  xlab("Position in the genome (Mb)") + 
  ylab("Coverage density") + 
  theme_bw() +
  ggtitle(label = "Fol4287 YP (mycelium-passaged) abd YC (Conidia-passaged) samples Genome coverage to 2025 Assembly (Ayhan et al.)",
          subtitle = "Median every 10kb per Scaffold, each sample normalized by its median, normalized to the WT") +
  theme(axis.text.x = element_text(size=6),
        strip.text.x = element_text(size = 8),
        axis.text.y = element_text(size=6)) +
  scale_x_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 1)) +
  scale_color_manual(values = c(
    "SRR7690004_cov" = "black",  
    "YC1_cov" = "#90BFF9", "YC2_cov" = "#90BFF9", "YC3_cov" = "#90BFF9", "YC4_cov" = "#90BFF9", "YC5_cov" = "#90BFF9",  
    "YP1_cov" = "#DE8BF9", "YP2_cov" = "#DE8BF9", "YP3_cov" = "#DE8BF9", "YP4_cov" = "#DE8BF9", "YP5_cov" = "#DE8BF9"   
  ))

ggsave("CNVS_median_MedianNorm_NormWT_pseudochromosomes.pdf", width = 10, height = 5, units = "in", dpi = 600)
