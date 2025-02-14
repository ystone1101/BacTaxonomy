# Load required libraries
library(ggplot2)
library(ggtree)
library(tidyr)
library(dplyr)
library(ggstar)  # Library for adding custom symbols to ggplot
library(ggnewscale)  # Allows multiple fill/colour scales in ggplot
library(ggtreeExtra)  # Extends ggtree for additional annotations
library(ggstance)  # Horizontal geoms for ggplot
library(patchwork)  # Library for combining multiple ggplot objects

# Load phylogenetic tree
data_dir <- "path/to/data"
plot_dir <- "path/to/output"

tree <- read.tree(file.path(data_dir, "phylogenomic_tree.nwk"))

# Generate tree plot
# The tree range is adjusted because multiple datasets will be attached on the left
p <- ggtree(tree) +
  geom_tiplab(align=TRUE, linesize=.5, linetype = NA, size=1) + xlim(0, 1)

# Load KEGG module data (Users can modify this dataset as needed)
# This dataset is based on METABOLIC (METABOLIC-G.pl)
genotype <- read.csv(file.path(data_dir, "keggmodule.csv"), stringsAsFactors=FALSE, check.names = FALSE)
genotype2 <- genotype %>% pivot_longer(-Module.ID, names_to = "variable", values_to = "value")

# Load antibiotic resistance genes (ARGs) data (Users can modify this dataset as needed)
args <- read.csv(file.path(data_dir, "ARGs.csv"), stringsAsFactor=F, check.names = F)
args2 <- args %>% pivot_longer(-ID, names_to = "variable", values_to = "value")

# Load peptidase enzyme data (Users can modify this dataset as needed)
peptidase <- read.csv(file.path(data_dir, "Peptidase2.csv"), stringsAsFactor=F, check.names = F)

# Load CAZyme (Carbohydrate-active enzymes) data (Users can modify this dataset as needed)
# This dataset is based on CAZy DB Hmmer results (Hidden Markov Models)
# There are 6 family groups in CAZy, but only 2 groups were detected in this study
cazy <- read.csv(file.path(data_dir, "CAZy.csv"), stringsAsFactor=F, check.names = F)

# Generate KEGG module heatmap
p1 <- ggplot(genotype2, aes(y=Module.ID, x=variable, fill=value)) + 
  geom_tile(color = "white",
            lwd = .5,
            linetype = 1) + 
  theme_classic() + 
  scale_x_discrete(position = "top") +
  scale_fill_manual(breaks=c("Present", "Absent"),  # Manually assign colors to categories 
                    values=c("lightblue", "grey90"), name="genotype") +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size =3, angle=90, hjust = -1, vjust = 1))
p1

# Generate antibiotic resistance genes heatmap
p2 <- ggplot(args2, aes(y=ID, x=variable, fill=value)) + 
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) + 
  theme_void() + # Removes all axes and background for a cleaner heatmap 
  scale_fill_manual(breaks=c("Present", "Absent"), 
                    values=c("lightgreen", "grey90"), name="genotype") +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())

p2

# Generate peptidase enzyme bar plot
p3 <- ggplot(peptidase, aes(ID, Gene)) + geom_col(aes(fill=Type), width = .75) + 
  coord_flip() + theme_tree2() + # Flips the X and Y axes, useful for horizontal bar plots 
  scale_y_continuous(limits = c(0,160), expand = c(0,1)) +
  scale_fill_manual(breaks=c("Aspartic", "Cycteine", "Glutamic", "Metallo", "Serine", "Threonine"),
                    values=c("#01B8AA", "#5F6B6D", "#A66999", "#FD625E", "#F2C80F", "#8AD4EB")) + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(color = "black", linewidth = .5))

p3

# Generate CAZyme enzyme bar plot
# Glycoside Hydrolases(GHs), GlycosylTransferases (GTs), Polysaccharide Lyases (PLs),  Auxiliary Activities (AAs), Carbohydrate Esterases (CEs), Carbohydrate-Binding Modules (CBMs)
# Full color set c("#0f2c40", "#1e3d5b", "#2f5b7e", "#4f7a8b", "#7fa0b1", "#a3c1e0")
p4 <- ggplot(cazy, aes(ID, Gene)) + 
  geom_col(aes(fill=Type), width = .75) + 
  coord_flip() + theme_tree2() + 
  scale_y_continuous(limits = c(0,300), expand = c(0,2)) +
  scale_fill_manual(breaks=c("GH", "PL"),  
                    values=c("#0f2c40", "#a3c1e0")) + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(color = "black", linewidth = .5))

p4

# Save KEGG module figure
# Inserts the tree plot to the left of the main plot
tiff(file.path(plot_dir, "Modules.tiff"), units="cm", width=16, height=7, res=900, compression = "lzw")
p5 <- p1 %>% insert_left(p,width = .3) 
p5
dev.off()

# Save characteristic genes figure
# Inserts the tree plot to the left of the main plot
# Adds space between plots for better readability
tiff(file.path(plot_dir, "charicteristic_genes.tiff"), units="cm", width=16, height=7, res=900, compression = "lzw")
p6 <- p2 %>% insert_left(p,width = .3) %>% insert_right(p3) %>% insert_right(plot_spacer(), width = .05) %>% insert_right(p4)
p6
dev.off()