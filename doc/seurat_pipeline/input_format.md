## Input file format
This section will describe all the parameters that can be set within the input parameters file. For each parameter, an example will be given in YAML format. A full template file can be found in [input_template.md](input_template.md).

#### Table of contents
- [project_name](#project_name)
- [output_path](#output_path)
- [species](#species)
- [data](#data)
- [include](#include)
- [nCells](#ncells)
- [nGenes](#ngenes)
- [mito_cutoff](#mito_cutoff)
- [batches](#batches)
- [variable_genes](#variable_genes)
- [pcs](#pcs)
- [point_size](#point_size)
- [tsne_plot](#tsne_plot)
- [feature_plot](#feature_plot)
- [resolution](#resolution)
- [violin_plot](#violin_plot)
- [dot_plot](#dot_plot)

##### project_name
This parameter is simply a string describing the dataset. It will only be used in the project.name slot of the Seurat object.
```yaml
project_name: my_project_name
```

##### output_path
This parameter is a path pointing to the folder where output files should be written to. No trailing / is needed.
```yaml
output_path: /path/to/output
```

##### species
This parameter should be set to either mouse or human depending on what species the samples are. This is used to identify mitochondrial genes, as mouse mitochondrial genes start with mt- and human mitochondrial genes start with MT-. This parameter is case insensitive and will cause an error if it does not match either mouse or human.
```yaml
species: mouse
```

##### data
This parameter contains maps listing all DGE files to be used with any relevant metadata. The name for each map will be used to set the orig.ident field in the Seurat object. Each map must contain a key path whose value is the path to the DGE file. Additional key value pairs can be specified which will be added to the metadata slot of the Seurat object, with the key indicating the metadata name and the value indicating the value of that metadata for the DGE. If a key appears for some DGEs but not all, it will be set to NA for all DGEs that lack the key. If the path key is not set or points to a non-existent file, the script will error.
```yaml
data:
  sample1:
    path: /path/to/sample1.dge.txt
    sex: Male
    diet: Fed
  sample2:
    path: /path/to/sample2.dge.txt
    diet: Fasted
```

##### include
This parameter is used to filter the dataset to only include certain cells. This should be a path to a file that contains one cell name per line. Multiple files can be given in a list. Any cell names that appear in at least one of the files will be kept. Cell names must match exactly. If a cell name in one of the files does not match any cells in the dataset, it will be ignored with no warning. This parameter is optional.
```yaml
include:
  - /path/to/good_cells_1.txt
  - /path/to/good_cells_2.txt
```

##### nCells
This parameter is used to set the minimum number of cells a gene must be present in. All genes that are present in less than this many cells are filtered out. To check multiple different cutoffs, set it to a list. If an element of the list is not numeric, the script will error.
```yaml
nCells:
  - 2
  - 5
```

##### nGenes
This parameter is used to set the minimum number of genes a cell must contain. All cells expressing fewer than this number of genes are filtered out. To check multiple different cutoffs, set it to a list. If an element of the list is not numeric, the script will error.
```yaml
nGenes:
  - 200
  - 500
```

##### mito_cutoff
This parameter is used to set the maximum fraction of mitochondrial genes a cell may contain. All cells with a higher mitochondrial gene ratio are filtered out. This number should be a fraction ranging from 0 (remove if any mitochondrial genes are present) to 1 (do not remove any cells). To check multiple different cutoffs, set it to a list. If an element of the list is not numeric, the script will error.
```yaml
mito_cutoff:
  - 0.1
  - 0.25
```

##### batches
This parameter is used to specify which batches to regress on. It is a map containing lists. The name of each list is used in the output naming scheme. Each element of the lists should be a piece of metadata to regress on. Each list is regressed on separately, items within lists are regressed on together. Possible metadata to regress on are percent.mito, nGene, nUMI, orig.ident, or any piece of metadata supplied with the data parameter. To skip regressing, supply a list with a name but no elements.
```yaml
batches:
  mito_and_ident:
    - percent.mito
    - orig.ident
  no_regression:
```

##### variable_genes
This parameter is used to select genes for PC generation. [M3Drop](https://github.com/tallulandrews/M3D) is used to calculate significance values for each gene. If this parameter is less than 1, it is interpreted as an FDR threshold. If it is greater than or equal to 1, it is interpreted as a number of genes to select sorted by FDR. To check multiple different genesets, set it to a list. If an element of the list is not numeric, the script will error.
```yaml
variable_genes:
  - 1000
  - .05
```

##### pcs
This parameter is used to select the PCs used for tSNE visualization and clustering. If this parameter is less than 1, it is interpreted as a cumulative proportion of variance threshold. The first n PCs will be selected such that the cumulative proportion of variation of those n PCs is greater than the threshold. If this number is greater than or equal to 1, it is interpreted as the number of PCs to select. Only 50 PCs are generated, so the number must be less than or equal to 50 and the proportion of variance is relative to the variance of the first 50 PCs. To check multiple different sets of PCs, set it to a list. If an element of the list is not numeric, the script will error.
```yaml
pcs:
  - 15
  - 0.6
```

##### point_size
This parameter controls the size of points in tSNE plots. This parameter is optional and will be set to 1 if it is not present. If it is not numeric, the script will error.
```yaml
point_size: 0.1
```

##### tsne_plot
This parameter specifies categorical metadata to plot on tSNE plots. orig.ident is always plotted by default. To plot multiple pieces of metadata, set it to a list. If a piece of metadata is not found, it will be omitted and a warning will be printed. This parameter is optional.
```yaml
tsne_plot:
  - sex
  - diet
```

##### feature_plot
This parameter specifies numerical metadata or genes to plot on tSNE plots. To plot multiple features, set it to a list. To plot multiple features on the same page, set an element of the list to a list. If a feature is not found, it will be omitted and a warning will be printed. This parameter is optional.
```yaml
feature_plot:
  - gene1
  -
    - gene2
    - gene3
    - gene4
    - gene5
```

##### resolution
This parameter specifies the clustering resolution. To check multiple different resolutions, set it to a list. If an element of the list is not numeric, the script will error.
```yaml
resolution:
  - 1.2
  - 1.5
```

##### violin_plot
This parameter specifies numerical metadata or genes to plot on violin plots. The x-axis of the violin plot are the clusters. percent.mito, nUMI, and nGene are plotted by default. To plot multiple features, set it to a list. To plot multiple features on the same page, set an element of the list to a list. If a feature is not found, it will be omitted and a warning will be printed. This parameter is optional.
```yaml
violin_plot:
  - gene1
  -
    - gene2
    - gene3
    - gene4
    - gene5
```

##### dot_plot
This parameter specifies genes to plot on dot plots. The y-axis of the dot plot are the clusters. To plot multiple features, set it to a list. To plot multiple features on the same page, set an element of the list to a list. If a feature is not found, it will be omitted and a warning will be printed. This parameter is optional.
```yaml
dot_plot:
  - gene1
  -
    - gene2
    - gene3
    - gene4
    - gene5
```
