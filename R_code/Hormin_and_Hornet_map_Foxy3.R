## Author = Ana Rodríguez López
## Date.of.creation = 12.06.24

# Generation of a map with locations of Hornet and Hormin copies

library(tidyverse)

# 1. Read file with Hornet and Hormin positions.
TE_positions <- read.table("Hormin-Hornet_positions_Foxy3.txt", sep = "\t", header = T)
Chromosome_lengths <- read.table("Chromosome_sizes_Foxy3.txt", sep = "\t", header = T)

# 2. Define chromosome levels for plotting in the correct order
chr_levels <- c("Unpositioned", "Chr15", "Chr14", "Chr13", "Chr12", "Chr11", 
               "Chr10", "Chr9", "Chr8", "Chr7", "Chr5", "Chr4", "Chr3", "Chr2", "Chr1")

# 3. Prepare chromosome data for plotting (BACKGROUND)
df_scaffolds <- Chromosome_lengths %>%
  mutate(
    Chromosome = factor(Chromosome, levels = chr_levels),
    Type = as.factor(Type),
    # Convert to kb
    Start_in_chr = Start_in_chr / 1000,
    End_in_chr = End_in_chr / 1000
  )

# 4. Prepare TE data for plotting (LINES)
df_TEs <- TE_positions %>%
  mutate(
    Chromosome = factor(Chromosome, levels = chr_levels),
    TE = as.factor(TE),
    # Convert to kb
    TE_position = TE_position / 1000
  )

# 5. Generate the chromosome plot with TE positions
plot_genome <- ggplot() +
  
  # Layer 1: Scaffold background
  geom_rect(data = df_scaffolds,
            aes(xmin = Start_in_chr, xmax = End_in_chr,
                ymin = as.numeric(Chromosome) - 0.3, 
                ymax = as.numeric(Chromosome) + 0.3,
                fill = Type),
            color = "black",    
            linewidth = 0.3) +  
  
  # Layer 2: TE lines
  geom_segment(data = df_TEs,
               aes(x = TE_position, xend = TE_position,
                   y = as.numeric(Chromosome) - 0.4, 
                   yend = as.numeric(Chromosome) + 0.4,
                   color = TE),
               linewidth = 0.7, alpha = 0.8) +
  
  # Colors
  scale_fill_manual(values = c("Core" = "#D3D3D3", "Fast-core" = "#D3D3D3","Accessory" = "#737373")) +
  scale_color_manual(values = c("Hormin" = "#f48849", "Hornet" = "#44bf70")) +
  
  # Apply chromosome names
  scale_y_continuous(breaks = 1:length(chr_levels),
                     labels = chr_levels) +
  
  # Adjust scale on X axis
  scale_x_continuous(breaks = seq(0, max(df_scaffolds$End_in_chr, na.rm = TRUE), by = 500)) +

  # Labels
  labs(x = "Position (kb)", 
       y = "Chromosome",
       fill = "Genomic compartment", 
       color = "TE") +
  
  # Other settings
  theme_light(base_size = 7) +
  theme(
    text = element_text(size = 7, face = "plain"), 
    axis.text.y = element_text(size = 7, face = "plain", color = "black"), 
    axis.text.x = element_text(size = 7, face = "plain", angle = 45, hjust = 1, color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )
print(plot_genome)

ggsave("positions_plot_Foxy3.bmp", plot = plot_genome, width = 5, height = 4, dpi = 320)
