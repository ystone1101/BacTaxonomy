# Load required libraries
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(ggplotify)
library(ggpubr)

# Define plot directory
plot_dir <- "C:/Users/yston/Desktop/plots"
individual_plot_dir <- file.path(plot_dir, "individual_plots")
start_time <- Sys.time()  # Start timing the script
if (!dir.exists(plot_dir)) dir.create(plot_dir)
if (!dir.exists(individual_plot_dir)) dir.create(individual_plot_dir)

# Define data directory
data_dir <- "C:/Users/yston/Desktop/CJ426/data"

# Function to load data and create heatmap
plot_heatmap <- function(file, color, tag, tag_position) {
  data <- read.csv(file, header = TRUE, row.names = 1, check.names = FALSE, fill = TRUE)
  color_palette <- colorRampPalette(brewer.pal(9, color))(100)
  p <- pheatmap(data,
                color = color_palette,
                border_color = "black",
                show_colnames = FALSE, 
                show_rownames = TRUE,
                treeheight_row = 15,
                treeheight_col = 15,
                width = 6,
                height = 2.5)
  as.ggplot(p) + labs(tag = tag) + theme(plot.tag.position = tag_position)
}

# Function to save individual plot
save_individual_plot <- function(plot, config_name, file) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%OS3")  # Include milliseconds for uniqueness
  outfile <- file.path(individual_plot_dir, paste0(config_name, "_", tools::file_path_sans_ext(basename(file)), "_", timestamp, ".tiff"))
  tiff(outfile, units = "in", width = 7.5, height = 5, res = 900, compression = "lzw")
  print(plot)
  dev.off()
}

# Define mapping for file names, colors, and tags
plot_config <- list(
  ANI = list(pattern = "ANI", color = "Reds", tag = "A. Ortho average nucleotide identity (OrthoANI)", tag_position = c(0.29, 1.06)),
  AAI = list(pattern = "AAI", color = "Blues", tag = "B. Average amino acid identity (AAI)", tag_position = c(0.23, 1.06)),
  POCP = list(pattern = "POCP", color = "Greens", tag = "C. The percentage of conserved proteins (POCP)", tag_position = c(0.3, 1.06)),
  DDH = list(pattern = "DDH", color = "Purples", tag = "D. DDH", tag_position = c(0.08, 1.04))
)

# Detect files and generate plots
files <- list.files(data_dir, pattern = "^(ANI|AAI|POCP|DDH).*\\.csv$", full.names = TRUE, ignore.case = TRUE)
print(files)  # Check detected files

plots <- list()
plot_names <- list()
log_file <- file.path(plot_dir, "error_log.txt")  # Log file for errors

for (config_name in names(plot_config)) {
  config <- plot_config[[config_name]]
  matched_files <- files[grepl(config$pattern, files)]
  if (length(matched_files) > 0) {
    for (file in matched_files) {
      plot_start_time <- Sys.time()
      plot <- tryCatch({
        plot_heatmap(file, config$color, config$tag, config$tag_position)
      }, error = function(e) {
        message(paste("Error processing file:", file, "-", e$message))
        file_size <- file.info(file)$size
        plot_execution_time <- Sys.time() - plot_start_time
        write(paste(Sys.time(), "-", file, "- Size:", file_size, "bytes - Execution Time:", round(plot_execution_time, 2), "seconds -", conditionMessage(e)), file = log_file, append = TRUE)
        NULL
      })
    }

if (!is.null(plot)) {
  plot_index <- length(plots) + 1
  plots[[plot_index]] <- plot
  plot_names[[plot_index]] <- paste0(config_name, "_", basename(file))
  save_individual_plot(plot, config_name, file)
  plot_execution_time <- Sys.time() - plot_start_time
  message(paste("Processed", file, "in", round(plot_execution_time, 2), "seconds"))
}
  }else {
  message(paste("No matching file found for", config_name))
}
}

# Display a warning if no matching files are found
if (length(plots) == 0) {
  message("No matching files found. No plots were generated.")
}

# Save the merged plot
outfile <- file.path(plot_dir, "Genome_comparison_value_merge_v4.tiff")
plot_height <- 3.5 * length(plots) + 0.5 * (length(plots) - 1)

if (length(plots) > 0) {
  tiff(outfile, units = "in", width = 7.5, height = plot_height, res = 900, compression = "lzw")
  
  layout_matrix <- cbind(rep(NA, length(plots) * 2 + 1))
  layout_matrix[seq(2, length(layout_matrix) - 1, by = 2)] <- seq_along(plots)
  heights <- rep(c(0.2, 1), length.out = length(layout_matrix))
  
  file_info <- file.info(files)
  ordered_files <- files[order(file_info$ctime)]  # Sort plots by file creation time
  plots <- plots[order(unlist(plot_names))]  # Sort plots by their names for consistent order  # Sort plots by their names for consistent order  # Sort plots by their names for consistent order
  
  g <- grid.arrange(grobs = plots, ncol = 1,
                    layout_matrix = layout_matrix,
                    heights = heights)
  execution_time <- Sys.time() - start_time
  print(paste("Execution Time:", execution_time))
  print(g)
  dev.off()  # Close the TIFF device after plotting
} else {
  message("No plots to arrange. Skipping grid.arrange().")
}
