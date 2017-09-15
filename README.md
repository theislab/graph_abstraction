# Graph abstraction relates subgroups of single cells


### Minimal examples with known ground truth

In [*minimal_examples*](minimal_examples), we study minimal example datasets
with known ground truth. In particular, a dataset that contains a tree-like
continuous manifold and clusters. We discuss in detail how we measure
robustness, how to reproduce figures of the paper and how competing algorithms
perform.

<img src="./sim_data/figures/aga_draw_graph_fr.png" height="200"><img src="./sim_data/figures/aga_graph.png" height="200">

### Hematopoietic lineage trees

#### Data from [Paul *et al.* (2015)](http://doi.org/10.1016/j.cell.2015.11.01)

<img src="./paul15/figures/draw_graph_fr.png" height="200"><img src="./paul15/figures/aga_graph.png" height="220">

For details, see [here](./paul15/).

#### Data from [Nestorowa *et al.* (2016)](http://doi.org/10.1182/blood-2016-05-716480)

<img src="./nestorowa16/figures/draw_graph_fr.png" height="200"><img src="./nestorowa16/figures/aga_graph.png" height="220">

For details, see [here](./nestorowa16).


### Lineage tree of whole adult animal

<img src="./figs_planaria/tsne.png" height="350"><img src="./figs_planaria/aga_graph.svg" height="450">

For details, see [here](./planaria).


### PBMC cells

From the literature.

<img src="./PBMCs_LineageTree.jpg" height="220">

All of the following datasets contain correct motifs from this tree, but as they are quite strongly clustering, the global tree cannot be inferred.

What would be needed would be data that contained all the intermediate states.

#### 3K cells from [10X (2017)](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k)

<img src="./figs_pbmc3k/tsne.png" height="200"><img src="./figs_pbmc3k/aga_graph.png" height="220">

For details, see [here](pbmc3k.html).

#### 33K cells from [10X (2017)](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc33k)

<img src="./figs_pbmc33k/tsne.png" height="200"><img src="./figs_pbmc33k/aga_graph.png" height="220">

For details, see [here](pbmc33k.html).

#### 68K cells from [Zheng *et al.* (2017)](https://doi.org/10.1038/ncomms14049)

<img src="./figs_zheng17/tsne.png" height="200"><img src="./figs_zheng17/aga_graph.png" height="220">

For details, see [here](zheng17.html).





