# Load required libraries
# Function to check and install missing packages
load_or_install <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}
load_or_install("ggplot2")
load_or_install("gridExtra")
load_or_install("grid")

# Define paths
data_path <- "/home/user/project/data/example_data.csv"
plot_dir <- "/home/user/project/plots"

# Create plot directory if it doesn't exist
plot_dir <- normalizePath(plot_dir, mustWork = FALSE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Load genome comparison data
Genome_comparison_parameter <- read.csv(data_path, header = TRUE, check.names = FALSE)
Genome <- as.data.frame(Genome_comparison_parameter)

# Define the %||% operator for default values
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Enhanced create scatter plot function
create_scatter_plot <- function(data, x_var, y_var, filename = NULL, hline = NA, vline = NA, annotations = NULL, tag = NULL, tag_position = NULL) {
  if (!(x_var %in% names(data))) stop(paste("Invalid x variable:", x_var))
  if (!(y_var %in% names(data))) stop(paste("Invalid y variable:", y_var))
  
  filename <- filename %||% paste0(x_var, "_vs_", y_var, ".tiff")
  tag_position <- tag_position %||% list(x = 0.2, y = 1.08)
  
  outfile <- file.path(plot_dir, filename)
  tiff(outfile, units="in", width=6, height=3, res=900, compression = "lzw")
  
  p <- ggplot(data, aes_string(x = x_var, y = y_var)) + 
    geom_point(size = 3.5, color = "black", aes(fill = Group, shape = Group)) +
    scale_shape_manual(values = c(22, 22)) +
    scale_fill_manual(values = c("gray", "red")) +
    guides(fill = FALSE, shape = FALSE) +
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"))
  
  # Conditional horizontal and vertical lines
  if (!is.na(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = "black")
  if (!is.na(vline)) p <- p + geom_vline(xintercept = vline, linetype = "dashed", color = "black")
  
  # Add custom annotations if provided
  if (!is.null(annotations) && !all(is.na(annotations))) {
    for (ann in annotations) {
      if (ann$geom == "text") {
        p <- p + annotate(geom = ann$geom, x = ann$x, y = ann$y,
                          label = ann$label, size = ann$size)
      } else if (ann$geom == "segment") {
        p <- p + annotate(geom = ann$geom, x = ann$x, y = ann$y,
                          xend = ann$xend, yend = ann$yend,
                          linewidth = ann$linewidth,
                          arrow = ann$arrow)
      }
    }
  }
  
  # Add tag if provided with customizable position
  if (!is.null(tag)) {
    p <- p + labs(tag = tag) + theme(plot.tag.position = c(tag_position$x, tag_position$y))
  }
  
  print(p)
  dev.off()
  return(p)
}

# Example annotations (modify as needed)
annotations1 <- list(
  list(geom = "text", x = 70.25, y = 88, label = "Inter genus", size = 3.5),
  list(geom = "text", x = 72.75, y = 88, label = "Intra genus", size = 3.5),
  list(geom = "segment", x = 72.5, xend = 69, y = 90, yend = 90, linewidth = 0.5, arrow = arrow(length = unit(0.25, 'cm'), type = 'closed')),
  list(geom = "segment", x = 72.5, xend = 74, y = 90, yend = 90, linewidth = 0.5, arrow = arrow(length = unit(0.25, 'cm'), type = 'closed'))
)

annotations2 <- list(
  list(geom = "text", x = 71.25, y = 78, label = "Different annotation", size = 3.5)
)

# Plot parameters
# Add or modify plot settings here to create more than 3 plots.
# Simply append additional list items with the desired plot configurations.
plot_params <- list(
  list(x_var = "orthoANI", y_var = "AAI", filename = NULL, 
       hline = NA, vline = 71.5, tag = "A. OrthoANI vs AAI", 
       tag_position = NULL, annotations = annotations1),
  
  list(x_var = "orthoANI", y_var = "POCP", filename = NULL, 
       hline = 60, vline = NA, tag = "B. OrthoANI vs POCP", 
       tag_position = NULL, annotations = annotations2),
  
  list(x_var = "AAI", y_var = "POCP", filename = NULL, 
       hline = NA, vline = NA, tag = "C. AAI vs POCP", 
       tag_position = NULL, annotations = NULL),
  
  list(x_var = "AAI", y_var = "POCP")
)

# Generate plots using lapply
plots <- lapply(plot_params, function(p) {
  create_scatter_plot(Genome, p$x_var, p$y_var, p$filename, 
                      hline = p$hline %||% NA, 
                      vline = p$vline %||% NA, 
                      annotations = p$annotations %||% NULL, 
                      tag = p$tag %||% NULL, 
                      tag_position = p$tag_position)
})

# Example Merged Plot
outfile <- file.path(plot_dir, "merge.tiff")
tiff(outfile, units = "in", width = 6, height = 12, res = 1200, compression = "lzw")

# Adjust figure spacing using layout_matrix and heights
layout_matrix <- cbind(rep(NA, length(plots) * 2 + 1))
layout_matrix[seq(2, length(layout_matrix) - 1, by = 2)] <- seq_along(plots)
heights <- rep(c(0.2, 1), length.out = nrow(layout_matrix))

g <- grid.arrange(grobs = plots, ncol = 1,
                  layout_matrix = layout_matrix,
                  heights = heights)

g
dev.off()