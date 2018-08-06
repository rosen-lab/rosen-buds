## Seurat Pipeline

This pipeline was designed for the [Rosen lab](http://www.evanrosenlab.net/) to facilitate the analysis of single cell RNAseq data using the [Seurat](http://satijalab.org/seurat/) R package. The R package [M3Drop](https://bioconductor.org/packages/release/bioc/html/M3Drop.html) is also used for finding variable genes. The pipeline goes from raw count matrices and outputs clustered and visualized data. It takes a single YAML file as input which specifies all the necessary parameters for the pipeline. Multiple parameters can be varied at once to test the effects of each parameter on the output.

The other markdown files in this directory provide instructions and information about running the pipeline:
- [running_the_pipeline.md](running_the_pipeline.md) goes through the pipeline script, explaining the inputs and outputs of each step.
- [input_format.md](input_format.md) describes the each parameter in the input YAML file.
- [input_template.md](input_template.md) contains a full template input YAML file.
- [yaml_primer.md](yaml_primer.md) gives a brief introduction to YAML syntax.
- [changelog.md](changelog.md) lists changes between different versions of the pipeline script.
