#! /broad/rosenlab_archive/Software/local/bin/Rscript

# Script to find and plot markers from Seurat pipeline output
version <- "0.2.1"

cat(format(Sys.time(), "[%H:%M:%S]"), "Loading libraries\n")

suppressMessages(library(methods))
suppressMessages(library(Seurat))
suppressMessages(library(yaml))
suppressMessages(library(viridis))

########################## FUNCTIONS ##########################
# Convenience function to timestamp output
cat_time <- function(x, ...) {
  cat(format(Sys.time(), "[%H:%M:%S]"), sprintf(x, ...))
}

# Function to load seurat object and rebuild its data slots
build_seurat <- function(path) {
  x <- readRDS(file.path(path, "seurat.rds"))
  x@data <- readRDS(file.path(path, "data.rds"))
  return(x)
}

# Return default parameters
get_default_params <- function() {
  return(
    list(
      ident.1 = "all",
      ident.2 = "all",
      logfc.threshold = 0.25,
      test.use = "wilcox",
      min.pct =  0.1,
      min.diff.pct = -Inf,
      only.pos = TRUE,
      return.thresh = 0.01,
      min.cells.gene = 3,
      min.cells.group = 3,
      latent.vars = "nUMI",
      name = "",
      pt.size = 1,
      tsne.pages = 1,
      tsne.genes.per.page = 9,
      dotplot.pages = 1,
      dotplot.genes.per.page = 9
    )
  )
}

# Function to load in parameters, replacing any missing parameters with defaults
load_params <- function(file) {
  # Get default parameters
  default <- get_default_params()

  # Return default params if file is not provided, otherwise read from yaml
  if (is.null(file)) {
    return(default)
  } else {
    params <- read_yaml(file)
  }

  # Fill in any missing parameters in yaml file
  for (param in names(default)) {
    if (is.null(params[[param]])) {
      params[[param]] <- default[[param]]
      cat(sprintf("Parameter %s not set, using default of %s\n", param, default[[param]]))
    }
  }

  return(params)
}

# Function to make and write tSNE plots
make_tSNE_plots <- function(split_markers, params, path) {
  pdf(file.path(path, sprintf("tSNE_plots_markers%s.pdf", params$name)))
  for (clust in split_markers) {
    for (page in 1:params$tsne.pages) {
      # Range of genes to plot on the page
      n_start <- (page-1)*params$tsne.genes.per.page + 1
      n_end <- n_start + params$tsne.genes.per.page - 1

      # Only plot if there are enough genes to plot
      n_end <- min(n_end, nrow(clust))
      if (n_start <= n_end) {
        FeaturePlot(seurat, clust$gene[n_start:n_end], no.legend = TRUE, cols.use = c("grey", "blue"), pt.size = params$pt.size)
      }
    }
  }
  invisible(dev.off())
}

# Function to make and write dot plots
make_dot_plots <- function(split_markers, params, path) {
  pdf(file.path(path, sprintf("dot_plots_markers%s.pdf", params$name)))
  TSNEPlot(seurat, do.label = TRUE, pt.size = params$pt.size)
  for (clust in split_markers) {
    for (page in 1:params$dotplot.pages) {
      # Range of genes to plot on the page
      n_start <- (page-1)*params$dotplot.genes.per.page + 1
      n_end <- n_start + params$dotplot.genes.per.page - 1

      # Only plot if there are enough genes to plot
      n_end <- min(n_end, nrow(clust))
      if (n_start <= n_end) {
        pdf(NULL)
        p <- DotPlot(seurat, clust$gene[n_start:n_end], do.return = TRUE, plot.legend = TRUE, x.lab.rot = TRUE)
        invisible(dev.off())
        suppressMessages(p <- p + scale_color_viridis(option = "plasma", direction = -1))
        p <- p + ggtitle(sprintf("Cluster %s: %i-%i", clust$cluster[1], n_start, n_end))
        print(p)
      }
    }
  }
  invisible(dev.off())
}

########################## SCRIPT ##########################
# Check command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat_time("No arguments provided, running with default parameters in current directory\n")
  path <- "."
  params_path <- NULL
} else if (length(args) == 1) {
  cat_time("Only one argument provided, running with default parameters\n")
  path <- args[1]
  params_path <- NULL
} else {
  cat_time("Loading parameters\n")
  path <- args[1]
  params_path <- args[2]
}

# Load params, add versions, and write params to file
params <- load_params(params_path)
params$script_version <- version
params$R_version <- R.version.string
params$Seurat_version <- paste(packageVersion("Seurat"), collapse = ".")
write_yaml(params, file.path(path, sprintf("markers%s.yml", params$name)))

cat_time("Loading Seurat object\n")
seurat <- build_seurat(path)

# Find markers
if (length(params$ident.2) == 1 && tolower(params$ident.2) == "all") {
  if (params$ident.1 == "all") {
    cat_time("Finding markers for all clusters\n")
    all_markers <- TRUE
  } else {
    cat_time("Finding markers for cluster %s against all clusters\n", params$ident.1)
    params$ident.2 <- NULL
    all_markers <- FALSE
  }
} else {
    if (length(params$ident.2) == 1) {
      cat_time("Finding markers for cluster %s against cluster %s\n", params$ident.1, params$ident.2)
      all_markers <- FALSE
    } else {
      cat_time("Finding markers for cluster %s against clusters %s\n", params$ident.1, paste(params$ident.2, collapse = ", "))
      all_markers <- FALSE
    }
}
if (all_markers) {
  markers <- FindAllMarkers(
    seurat,
    do.print = TRUE,
    logfc.threshold = params$logfc.threshold,
    test.use = params$test.use,
    min.pct = params$min.pct,
    min.diff.pct = params$min.diff.pct,
    only.pos = params$only.pos,
    return.thresh = params$return.thresh,
    min.cells.gene = params$min.cells.gene,
    min.cells.group = params$min.cells.group,
    latent.vars = params$latent.vars
  )
} else {
  markers <- FindMarkers(
    seurat,
    do.print = TRUE,
    ident.1 = params$ident.1,
    ident.2 = params$ident.2,
    logfc.threshold = params$logfc.threshold,
    test.use = params$test.use,
    min.pct = params$min.pct,
    min.diff.pct = params$min.diff.pct,
    only.pos = params$only.pos,
    return.thresh = params$return.thresh,
    min.cells.gene = params$min.cells.gene,
    min.cells.group = params$min.cells.group,
    latent.vars = params$latent.vars
  )
  markers$gene <- row.names(markers)
}

cat_time("Writing markers\n")
write.table(
  markers, file.path(path, sprintf("markers%s.tsv", params$name)),
  quote = FALSE, sep = "\t", row.names = FALSE
)

# Split markers by cluster, need to add cluster if FindMarkers was run
if (!all_markers) {
  markers$cluster <- as.factor(params$ident.1)
}
split_markers <- split(markers, markers$cluster)

cat_time("Making tSNE plots\n")
make_tSNE_plots(split_markers, params, path)

cat_time("Making dot plots\n")
make_dot_plots(split_markers, params, path)

cat_time("Done\n")
