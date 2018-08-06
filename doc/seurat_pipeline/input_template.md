## YAML Example
Below is a complete input template file. The file is in YAML syntax. A brief primer on YAML syntax can be found in [yaml_primer.md](yaml_primer.md). A description of each input parameter is can be found in [input_format.md](input_format.md).

```yaml
project_name: my_project_name

output_path: /path/to/output

species: mouse

data:
  sample1:
    path: /path/to/sample1.dge.txt
    sex: Male
    diet: Fed
  sample2:
    path: /path/to/sample2.dge.txt
    diet: Fasted
    
include:
  - /path/to/good_cells_1.txt
  - /path/to/good_cells_2.txt

nCells:
  - 2
  - 5

nGenes:
  - 200
  - 500

mito_cutoff:
  - 0.1
  - 0.25

batches:
  mito_and_ident:
    - percent.mito
    - orig.ident
  no_regression:

variable_genes:
  - 1000
  - .05

pcs:
  - 15
  - 0.6

point_size: 0.1

tsne_plot:
  - sex
  - diet

feature_plot:
  - gene1
  -
    - gene2
    - gene3
    - gene4
    - gene5

resolution:
  - 1.2
  - 1.5

violin_plot:
  - gene1
  -
    - gene2
    - gene3
    - gene4
    - gene5

dot_plot:
  - gene1
  -
    - gene2
    - gene3
    - gene4
    - gene5
```
