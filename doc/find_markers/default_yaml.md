## Default options
This page contains a YAML file with the default parameters for the Find Markers script. These parameters are used if no YAML file is provided or if individual parameters are missing.

```yaml
ident.1: all
ident.2: all

logfc.threshold: 0.25
test.use: wilcox
min.pct:  0.1
min.diff.pct: -Inf
only.pos: TRUE
return.thresh: 0.01
min.cells.gene: 3
min.cells.group: 3
latent.vars: ["nUMI"]

name: ""
pt.size: 1
tsne.pages: 1
tsne.genes.per.page: 9
dotplot.pages: 1
dotplot.genes.per.page: 9 

```