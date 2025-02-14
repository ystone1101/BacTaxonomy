library(plotrix)
library(wesanderson)

# Function to create a flower plot with adjustable transparency
flower_plot <- function(sample, value, start = 90, a = 2.5, b = 1.5,  
                        ellipse_col = rgb(135, 206, 235, 150, max = 255), 
                        circle_col = rgb(0, 162, 214, max = 255),
                        circle_text_cex = 1, labels = NULL, shared_gene_count = NULL) {
  
  # Function to adjust color transparency
  addalpha <- function(colors, alpha=1.0) {
    r <- col2rgb(colors, alpha=T)
    # Apply alpha
    r[4,] <- alpha * 255
    r <- r / 255.0
    return(rgb(r[1,], r[2,], r[3,], r[4,]))
  }
  
  # Set default label if NULL, including shared_gene_count
  if (is.null(labels)) {
    labels <- paste(" Core gene \n families", if (!is.null(shared_gene_count)) paste0("\n", shared_gene_count) else "")
  }
  
  # Ensure sample and value vectors match in length
  if (length(sample) != length(value)) {
    stop("Error: 'sample' and 'value' must have the same length.")
  }
  
  # Set up plot parameters (removes axis and margins)
  par(bty = "n", ann = FALSE, xaxt = "n", yaxt = "n", mar = c(0, 0, 0, 0))
  plot(c(0, 10), c(0, 10), type = "n")  # Empty plot to set coordinate space
  
  # Number of petals in the flower
  n <- length(sample)
  deg <- 360 / n  # Angle interval between petals
  
  # Ensure ellipse_col has enough colors
  ellipse_col <- rep(ellipse_col, length.out = n)
  
  # Draw petals
  lapply(1:n, function(t) {
    plotrix::draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180), 
                          y = 5 + sin((start + deg * (t - 1)) * pi / 180), 
                          col = ellipse_col[t],  # 색상 적용
                          border = "black",
                          a = a, b = b, angle = deg * (t - 1))
    
    # Display gene counts inside petals (Adjustable)
    text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
         y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
         value[t], font = 1)
    
    # Display sample labels around petals (Adjustable)
    rotation <- ifelse(deg * (t - 1) < 180, deg * (t - 1) - start, deg * (t - 1) + start)
    text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
         y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
         sample[t], srt = rotation, adj = ifelse(rotation < 180, 1, 0), cex = circle_text_cex)
  })
  
  # Draw center circle to represent core gene families
  plotrix::draw.circle(x = 5, y = 5, r = 1, col = "white", border = "black")
  
  # Display central label
  text(x = 5, y = 5, labels = labels, font = 1.5)  # Core gene families label
}


# Load gene presence/absence data
df <- read.csv("path/to/gf_data.csv", header = TRUE, check.names = FALSE)

if ("Gene" %in% colnames(df)) {
  df <- df[, !colnames(df) %in% "Gene"]
}  # Remove 'Gene' column if necessary

calculate_gene_counts <- function(df) {
  sample_names <- colnames(df)
  
  # Unique gene count (genes present in only one sample)
  unique_counts <- sapply(sample_names, function(sample) {
    other_samples <- setdiff(sample_names, sample)
    unique_genes <- df[[sample]] == 1 & rowSums(df[other_samples], na.rm = TRUE) == 0
    sum(unique_genes, na.rm = TRUE)
  })
  
  # Shared gene count (genes present in ALL samples)
  shared_genes <- rowSums(df, na.rm = TRUE) == ncol(df)
  shared_gene_count <- sum(shared_genes, na.rm = TRUE)
  
  # Return as separate data frames
  unique_gene_counts <- data.frame(Sample = seq_along(sample_names), Unique_Gene_Counts = unique_counts)
  shared_gene_counts <- data.frame(Shared_Gene_Count = shared_gene_count)
  
  return(list(unique_gene_counts = unique_gene_counts, shared_gene_counts = shared_gene_counts))
}

# Calculate the number of genes from the dataframe
gene_counts <- calculate_gene_counts(df)

# Extract unique and shared gene counts separately
unique_gene_counts <- gene_counts$unique_gene_counts
shared_gene_count <- gene_counts$shared_gene_counts$Shared_Gene_Count

# Display results
print(unique_gene_counts)
print(shared_gene_count)

if (length(unique_gene_counts$Sample) > length(LETTERS)) {
  warning("Too many samples! Consider splitting into multiple plots.")
}


tiff(
  filename="path/to/Venn_test.tiff",
  width=15,
  height=15,
  unit="cm",
  res = 900,
  compression = "lzw"
  )

par(mar = c(0, 0, 0, 0))

flower_plot(unique_gene_counts$Sample, 
            unique_gene_counts$Unique_Gene_Counts, 
            90, 
            0.8, 
            2,
            shared_gene_count = shared_gene_count,
            ellipse_col = addalpha(colorRampPalette(wes_palette("Darjeeling1"))(length(unique_gene_counts$Sample)), 0.5))

dev.off()