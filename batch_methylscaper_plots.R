#!/usr/bin/env Rscript
#' batch_methylscaper_plots.R
#' 
#' Batch generate methylscaper plots from bsFreqs CSV output files
#' 
#' Usage:
#'   Rscript batch_methylscaper_plots.R <input_dir> [output_dir] [--width 8] [--height 10] [--format png]
#'
#' Arguments:
#'   input_dir     Directory containing *_map.csv files from bsFreqs.py --csv
#'   output_dir    Directory for output plots (default: input_dir/plots)
#'   --left-site   Site for left panel (default: auto-detect HCG/CG)
#'   --right-site  Site for right panel (default: auto-detect GCH/GC)
#'   --width       Plot width in inches (default: 8)
#'   --height      Plot height in inches (default: 10)
#'   --format      Output format: png, pdf, or both (default: png)
#'   --dpi         Resolution for PNG output (default: 150)
#'
#' Example:
#'   Rscript batch_methylscaper_plots.R /path/to/frequencies/
#'   Rscript batch_methylscaper_plots.R /path/to/frequencies/ --left-site CG --right-site GC
#'
#' Author: Jason Orr Brant, 2026

# ==============================================================================
# Load required functions
# ==============================================================================

#' Load bsFreqs CSV files directly into plotting format
#' (Embedded from load_bsFreqs_for_methylscaper.R)
load_bsFreqs_csv_direct <- function(hcg_file, gch_file, method = "PCA") {
  
  # Read files - keeping first column as row names
  hcg_raw <- read.table(hcg_file, 
                        header = TRUE, 
                        sep = "\t", 
                        comment.char = "",
                        check.names = FALSE,
                        row.names = 1,
                        stringsAsFactors = FALSE,
                        na.strings = ".")
  
  gch_raw <- read.table(gch_file, 
                        header = TRUE, 
                        sep = "\t", 
                        comment.char = "",
                        check.names = FALSE,
                        row.names = 1,
                        stringsAsFactors = FALSE,
                        na.strings = ".")
  
  # Convert to numeric matrices
  hcg_matrix <- apply(hcg_raw, 2, as.numeric)
  gch_matrix <- apply(gch_raw, 2, as.numeric)
  
  # Handle case where apply returns a vector (single row)
  if (is.null(dim(hcg_matrix))) {
    hcg_matrix <- matrix(hcg_matrix, nrow = 1)
    gch_matrix <- matrix(gch_matrix, nrow = 1)
  }
  
  # Replace NA with 0 (will show as white/background after recode)
  hcg_matrix[is.na(hcg_matrix)] <- 0
  gch_matrix[is.na(gch_matrix)] <- 0
  
  # Apply the recode transformation
  # GCH transformation (for yellow/black display on RIGHT side of plot)
  gch_recoded <- gch_matrix
  gch_recoded[gch_matrix == 2] <- -4
  gch_recoded[gch_matrix == 1] <- -3
  gch_recoded[gch_matrix == 0] <- -2.5
  gch_temp_neg1 <- gch_matrix == -1
  gch_recoded[gch_matrix == -2] <- -1
  gch_recoded[gch_temp_neg1] <- -2
  gch_recoded[gch_matrix == 3] <- 5 # maps to blue for wrong bases
  gch_recoded[gch_matrix == 4] <- 4.5 # fill between wrong bases
  
  # HCG transformation (for red/black display on LEFT side of plot)
  hcg_recoded <- hcg_matrix
  hcg_recoded[hcg_matrix == 2] <- 4
  hcg_recoded[hcg_matrix == 1] <- 3
  hcg_recoded[hcg_matrix == 0] <- 2.5
  hcg_recoded[hcg_matrix == -1] <- 2
  hcg_recoded[hcg_matrix == -2] <- 1
  hcg_recoded[hcg_matrix == 3] <- 5
  hcg_recoded[hcg_matrix == 4] <- 4.5
  
  # Combine into toClust format: GCH columns first, then HCG columns
  toClust <- cbind(gch_recoded, hcg_recoded)
  
  # Ordering by first principal component
  if (method == "PCA" && nrow(toClust) > 1) {
    col_centered <- apply(toClust, 2, function(x) x - mean(x, na.rm = TRUE))
    col_centered[is.na(col_centered)] <- 0
    svd_result <- svd(col_centered, nu = 1, nv = 0)
    order1 <- order(svd_result$u[, 1])
  } else {
    order1 <- seq_len(nrow(toClust))
  }
  
  return(list(
    toClust = toClust,
    order1 = order1
  ))
}


#' Generate methylation sequence plot
#' (Adapted from methylscaper's plotSequence function)
plotSequence <- function(orderObject, plotFast = TRUE,
                         blankWidth = NULL, Title = "",
                         drawLine = TRUE, drawKey = TRUE,
                         leftLabel = "HCG", rightLabel = "GCH") {
  
  toClust <- orderObject$toClust
  order1 <- orderObject$order1
  input_GCH <- toClust[, seq(1, (ncol(toClust) / 2)), drop = FALSE]
  input_HCG <- toClust[, seq((ncol(toClust) / 2 + 1), ncol(toClust)), drop = FALSE]
  
  myCols <- c(
    "darkgoldenrod2", "yellow", "gray62", "black",
    "gray80", "white", "gray80", "black", "gray62", 
    "red", "darkred", "lightblue", "blue", "blue"
  )
  myBreaks <- c(-5, -4, -3, -2.5, -2, -1, 0, 1, 2, 2.5, 3, 4, 4.5, 5, 6)
  
  if (nrow(toClust) == 1) {
    input_HCG <- t(as.matrix(input_HCG))
    input_GCH <- t(as.matrix(input_GCH))
  }
  
  input_HCG_fix <- input_HCG
  input_GCH_fix <- input_GCH
  
  if (is.null(blankWidth)) blankWidth <- round(.12 * ncol(toClust) / 2)
  
  blankCOLS <- matrix(rep(0, nrow(input_HCG) * blankWidth),
                      nrow = nrow(input_HCG), ncol = blankWidth)
  toPlot_fix_og <- cbind(input_HCG_fix, blankCOLS, input_GCH_fix)
  
  sites <- which(apply(abs(toPlot_fix_og), 2, function(x) any(x %in% c(4, 1))))
  sites_scale <- sites / ncol(toPlot_fix_og)
  
  toPlot_fix <- toPlot_fix_og[rev(order1), , drop = FALSE]
  
  if (nrow(toClust) == 1) {
    toPlot_fix <- t(as.matrix(toPlot_fix))
  }
  
  # Add legend key
  if (drawKey == TRUE) {
    totalHeight <- pmax(ceiling(nrow(toPlot_fix) * 0.1), 3)
    blankROW <- matrix(rep(0, ncol(toPlot_fix) * totalHeight),
                       nrow = totalHeight, ncol = ncol(toPlot_fix))
    keyHeight <- pmax(ceiling(totalHeight * 0.5), 2)
    blankROW[seq(1, keyHeight), seq(1, min(147, ncol(input_HCG_fix)))] <- 2
    blankROW[seq(1, keyHeight), (seq(1, min(147, ncol(input_GCH_fix))) + 
                                   (ncol(input_HCG_fix) + ncol(blankCOLS)))] <- 2
    toPlot_fix <- rbind(blankROW, toPlot_fix)
  }
  
  # Plotting
  par(xpd = FALSE, mar = c(2, 2, 4, 1), mgp = c(0, 0.5, 0))
  image(t(toPlot_fix),
        col = myCols, axes = FALSE, breaks = myBreaks,
        useRaster = plotFast, ylim = c(0, 1))
  
  # Site tick marks
  axis(3,
       at = sites_scale, labels = rep("", length(sites_scale)),
       tick = TRUE, line = 0.5, col = "white", cex = 1, lwd = 1,
       col.ticks = "black", tck = -0.02)
  
  title(Title, line = 2, col.main = "royalblue4")
  
  # Axis labels
  plot1 <- round(ncol(input_HCG_fix) / ncol(toPlot_fix), 2)
  plot2 <- round((ncol(input_HCG_fix) + blankWidth) / ncol(toPlot_fix), 2)
  
  title(leftLabel, adj = plot1 / 2 - 0.025, line = 1.5)
  title(rightLabel, adj = plot2 + 0.025 + (1 - plot2) / 2, line = 1.5)
  
  toLabel <- rev(seq(1, length(order1), by = ceiling(length(order1) / 8)))
  if (!(length(order1) %in% toLabel)) toLabel <- c(length(order1), toLabel)
  y_axis_starting_point <- ifelse(drawKey, nrow(blankROW) / nrow(toPlot_fix), 0)
  axis(2, at = seq(y_axis_starting_point, 1, length.out = length(toLabel)), 
       labels = toLabel)
  
  toLabel <- round(c(
    seq(1, ncol(input_HCG_fix), length.out = 5),
    seq(1, ncol(input_GCH_fix), length.out = 5)
  ))
  axis(1, at = seq(0, plot1, length.out = 5), labels = toLabel[seq(1, 5)])
  axis(1, at = seq(plot2, 1, length.out = 5), labels = toLabel[seq(1, 5)])
  
  # DNA line
  if (drawLine == TRUE) {
    par(xpd = NA)
    top1 <- 1.01
    segments(0, top1, plot1, top1, lwd = 1)
    segments(plot2, top1, 1, top1, lwd = 1)
  }
}


# ==============================================================================
# Main batch processing function
# ==============================================================================

#' Find matching site file pairs in a directory
#' 
#' Looks for paired CSV map files for the left-panel site (endogenous/red)
#' and right-panel site (accessibility/yellow).
#' 
#' @param input_dir Directory containing *_map.csv files
#' @param left_site Site name for left panel (default: auto-detect)
#' @param right_site Site name for right panel (default: auto-detect)
#' @return List of file pairs with hcg, gch, identifier, left_site, right_site
find_csv_pairs <- function(input_dir, left_site = NULL, right_site = NULL) {
  
  # Find all *_map.csv files, excluding *_coded_map.csv (combined format)
  all_map_files <- list.files(input_dir, pattern = "_map\\.csv$", full.names = TRUE)
  all_map_files <- all_map_files[!grepl("_coded_map\\.csv$", all_map_files)]
  
  if (length(all_map_files) == 0) {
    return(list())
  }
  
  # Extract site types from filenames: {SITE}-{strand}-{locus}_map.csv
  file_info <- lapply(all_map_files, function(f) {
    bn <- basename(f)
    # Remove _map.csv suffix
    stem <- sub("_map\\.csv$", "", bn)
    # Split on first dash to get site type
    parts <- strsplit(stem, "-", fixed = TRUE)[[1]]
    if (length(parts) >= 2) {
      site <- parts[1]
      identifier <- paste(parts[-1], collapse = "-")
      return(list(file = f, site = site, identifier = identifier))
    }
    return(NULL)
  })
  file_info <- Filter(Negate(is.null), file_info)
  
  if (length(file_info) == 0) {
    return(list())
  }
  
  # Get unique site types present
  available_sites <- unique(sapply(file_info, function(x) x$site))
  cat("  Available site types in files:", paste(available_sites, collapse = ", "), "\n")
  
  # Auto-detect site roles if not specified
  # Left panel = endogenous (red): HCG, CG, WCG, or first site alphabetically
  # Right panel = accessibility (yellow): GCH, GC, or second site
  endogenous_sites <- c("HCG", "CG", "WCG", "SCG")  # Sites measuring endogenous methylation
  accessibility_sites <- c("GCH", "GC")               # Sites measuring accessibility
  
  if (is.null(left_site)) {
    # Try to find an endogenous site
    found_endo <- intersect(available_sites, endogenous_sites)
    if (length(found_endo) > 0) {
      left_site <- found_endo[1]
    } else if (length(available_sites) >= 2) {
      left_site <- sort(available_sites)[1]
    } else {
      left_site <- available_sites[1]
    }
  }
  
  if (is.null(right_site)) {
    # Try to find an accessibility site
    found_acc <- intersect(available_sites, accessibility_sites)
    if (length(found_acc) > 0) {
      right_site <- found_acc[1]
    } else if (length(available_sites) >= 2) {
      remaining <- setdiff(available_sites, left_site)
      right_site <- sort(remaining)[1]
    } else {
      right_site <- NULL
    }
  }
  
  cat("  Left panel (endogenous/red):", left_site, "\n")
  cat("  Right panel (accessibility/yellow):", ifelse(is.null(right_site), "NONE", right_site), "\n")
  
  # If only one site type exists, we can't make paired plots
  if (is.null(right_site) || !(left_site %in% available_sites) || !(right_site %in% available_sites)) {
    cat("  WARNING: Need exactly two site types for paired plots.\n")
    cat("  Available:", paste(available_sites, collapse = ", "), "\n")
    cat("  Use --left-site and --right-site to specify which sites to pair.\n")
    return(list())
  }
  
  # Build lookup: identifier -> file, by site type
  left_files <- list()
  right_files <- list()
  
  for (fi in file_info) {
    if (fi$site == left_site) {
      left_files[[fi$identifier]] <- fi$file
    } else if (fi$site == right_site) {
      right_files[[fi$identifier]] <- fi$file
    }
  }
  
  # Match pairs by identifier
  pairs <- list()
  common_ids <- intersect(names(left_files), names(right_files))
  
  for (id in common_ids) {
    pairs[[id]] <- list(
      hcg = left_files[[id]],      # "hcg" = left panel (endogenous)
      gch = right_files[[id]],      # "gch" = right panel (accessibility)
      identifier = id,
      left_site = left_site,
      right_site = right_site
    )
  }
  
  # Warn about unmatched files
  unmatched_left <- setdiff(names(left_files), common_ids)
  unmatched_right <- setdiff(names(right_files), common_ids)
  if (length(unmatched_left) > 0) {
    for (id in unmatched_left) {
      warning(paste("No matching", right_site, "file found for:", left_site, id))
    }
  }
  if (length(unmatched_right) > 0) {
    for (id in unmatched_right) {
      warning(paste("No matching", left_site, "file found for:", right_site, id))
    }
  }
  
  return(pairs)
}


#' Process a single file pair and generate plot
process_pair <- function(pair, output_dir, width, height, format, dpi) {
  
  identifier <- pair$identifier
  left_site <- if (!is.null(pair$left_site)) pair$left_site else "HCG"
  right_site <- if (!is.null(pair$right_site)) pair$right_site else "GCH"
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
  cat("Processing:", identifier, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")
  
  cat("  Left panel (", left_site, "):", basename(pair$hcg), "\n", sep = "")
  cat("  Right panel (", right_site, "):", basename(pair$gch), "\n", sep = "")
  
  # Load data
  tryCatch({
    orderObj <- load_bsFreqs_csv_direct(pair$hcg, pair$gch, method = "PCA")
    
    cat("  Sequences:", nrow(orderObj$toClust), "\n")
    cat("  Positions:", ncol(orderObj$toClust) / 2, "per panel\n")
    
    # Generate plot title from identifier
    # identifier format: {strand}-{locus}
    parts <- strsplit(identifier, "-", fixed = TRUE)[[1]]
    strand <- parts[1]
    locus <- paste(parts[-1], collapse = "-")
    plot_title <- paste0(locus, " (", ifelse(strand == "A", "Top", "Bottom"), " strand)")
    
    # Output filename base
    out_base <- file.path(output_dir, paste0("methylmap_", identifier))
    
    # Generate PNG
    if (format %in% c("png", "both")) {
      png_file <- paste0(out_base, ".png")
      png(png_file, width = width, height = height, units = "in", res = dpi)
      plotSequence(orderObj, Title = plot_title, plotFast = TRUE,
                   leftLabel = left_site, rightLabel = right_site)
      dev.off()
      cat("  Wrote:", basename(png_file), "\n")
    }
    
    # Generate PDF
    if (format %in% c("pdf", "both")) {
      pdf_file <- paste0(out_base, ".pdf")
      pdf(pdf_file, width = width, height = height)
      plotSequence(orderObj, Title = plot_title, plotFast = FALSE,
                   leftLabel = left_site, rightLabel = right_site)
      dev.off()
      cat("  Wrote:", basename(pdf_file), "\n")
    }
    
    return(TRUE)
    
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    return(FALSE)
  })
}


#' Main batch processing function
batch_process <- function(input_dir, output_dir = NULL, 
                          width = 8, height = 10, 
                          format = "png", dpi = 150,
                          left_site = NULL, right_site = NULL) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║     Methylscaper Batch Plot Generator                    ║\n")
  cat("║     Kladde Lab, 2026                                     ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Validate input directory
  if (!dir.exists(input_dir)) {
    stop(paste("Input directory does not exist:", input_dir))
  }
  
  # Set default output directory
  if (is.null(output_dir)) {
    output_dir <- file.path(input_dir, "plots")
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  cat("Input directory:", input_dir, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Plot dimensions:", width, "x", height, "inches\n")
  cat("Output format:", format, "\n")
  if (format %in% c("png", "both")) {
    cat("PNG resolution:", dpi, "dpi\n")
  }
  
  # Find CSV pairs
  pairs <- find_csv_pairs(input_dir, left_site = left_site, right_site = right_site)
  
  if (length(pairs) == 0) {
    cat("\nNo matching site file pairs found in:", input_dir, "\n")
    cat("Expected file pattern: {SITE}-{strand}-{locus}_map.csv\n")
    cat("Example: HCG-A-locus_map.csv and GCH-A-locus_map.csv\n")
    cat("         or: CG-A-locus_map.csv and GC-A-locus_map.csv\n")
    cat("\nUse --left-site and --right-site to specify which sites to pair.\n")
    return(invisible(FALSE))
  }
  
  cat("\nFound", length(pairs), "file pair(s) to process\n")
  
  # Process each pair
  success_count <- 0
  fail_count <- 0
  
  for (pair in pairs) {
    result <- process_pair(pair, output_dir, width, height, format, dpi)
    if (result) {
      success_count <- success_count + 1
    } else {
      fail_count <- fail_count + 1
    }
  }
  
  # Summary
  cat("\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("SUMMARY\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("  Successful:", success_count, "\n")
  cat("  Failed:", fail_count, "\n")
  cat("  Output directory:", output_dir, "\n")
  cat("\n")
  
  return(invisible(success_count > 0))
}


# ==============================================================================
# Command line interface
# ==============================================================================

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Defaults
  params <- list(
    input_dir = NULL,
    output_dir = NULL,
    width = 8,
    height = 10,
    format = "png",
    dpi = 150,
    left_site = NULL,
    right_site = NULL
  )
  
  # Show help if no arguments
  if (length(args) == 0 || args[1] %in% c("-h", "--help")) {
    cat("
Usage: Rscript batch_methylscaper_plots.R <input_dir> [output_dir] [options]

Arguments:
  input_dir       Directory containing *_map.csv files from bsFreqs.py --csv
  output_dir      Directory for output plots (default: input_dir/plots)

Options:
  --left-site S   Site type for left panel / endogenous (default: auto-detect HCG or CG)
  --right-site S  Site type for right panel / accessibility (default: auto-detect GCH or GC)
  --width N       Plot width in inches (default: 8)
  --height N      Plot height in inches (default: 10)
  --format F      Output format: png, pdf, or both (default: png)
  --dpi N         Resolution for PNG output (default: 150)
  -h, --help      Show this help message

Site auto-detection:
  The script automatically detects which site types are available in the
  input directory and assigns them to left/right panels:
    Left panel (endogenous/red):     HCG > CG > WCG > first alphabetically
    Right panel (accessibility/yellow): GCH > GC > second alphabetically

Examples:
  # Standard HCG/GCH analysis (auto-detected)
  Rscript batch_methylscaper_plots.R /path/to/frequencies/

  # Explicit CG/GC site pairing
  Rscript batch_methylscaper_plots.R /path/to/frequencies/ --left-site CG --right-site GC

  # Custom output directory and format
  Rscript batch_methylscaper_plots.R /path/to/frequencies/ /path/to/plots --format both --dpi 300

")
    quit(status = 0)
  }
  
  # Parse positional and optional arguments
  i <- 1
  positional <- c()
  
  while (i <= length(args)) {
    arg <- args[i]
    
    if (arg == "--left-site" && i < length(args)) {
      params$left_site <- args[i + 1]
      i <- i + 2
    } else if (arg == "--right-site" && i < length(args)) {
      params$right_site <- args[i + 1]
      i <- i + 2
    } else if (arg == "--width" && i < length(args)) {
      params$width <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (arg == "--height" && i < length(args)) {
      params$height <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (arg == "--format" && i < length(args)) {
      params$format <- args[i + 1]
      i <- i + 2
    } else if (arg == "--dpi" && i < length(args)) {
      params$dpi <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (!startsWith(arg, "-")) {
      positional <- c(positional, arg)
      i <- i + 1
    } else {
      cat("Unknown option:", arg, "\n")
      i <- i + 1
    }
  }
  
  # Assign positional arguments
  if (length(positional) >= 1) {
    params$input_dir <- positional[1]
  }
  if (length(positional) >= 2) {
    params$output_dir <- positional[2]
  }
  
  # Validate
  if (is.null(params$input_dir)) {
    stop("Input directory is required. Use -h for help.")
  }
  
  if (!params$format %in% c("png", "pdf", "both")) {
    stop("Format must be 'png', 'pdf', or 'both'")
  }
  
  return(params)
}


# ==============================================================================
# Main entry point
# ==============================================================================

main <- function() {
  params <- parse_args()
  
  success <- batch_process(
    input_dir = params$input_dir,
    output_dir = params$output_dir,
    width = params$width,
    height = params$height,
    format = params$format,
    dpi = params$dpi,
    left_site = params$left_site,
    right_site = params$right_site
  )
  
  quit(status = if (success) 0 else 1)
}

# Run if called as script
if (!interactive()) {
  main()
}
