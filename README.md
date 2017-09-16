# Graph abstraction relates subgroups of single cells


### Minimal examples with known ground truth

In [*minimal_examples*](minimal_examples), we study clean simulated datasets
with known ground truth. In particular, a dataset that contains a tree-like
continuous manifold and clusters. We discuss in detail how we measure
robustness, how to reproduce figures of the paper and how competing algorithms
perform.

<img src="./minimal_examples/figures/aga_draw_graph_fr.png" height="200"><img src="./minimal_examples/figures/aga_graph.png" height="200">

### Hematopoietic lineage trees

Here, we consider two well-studied datasets on hematopoietic differentiation.

#### Data from [Paul *et al.* (2015)](http://doi.org/10.1016/j.cell.2015.11.01)

In [*paul15*](paul15), we analyze data for myeloid progenitor development. This is the same
data, which has served as benchmark for Monocle 2 [(Qiu *et al.*,
  Nat. Meth., 2017)](https://doi.org/10.1038/nmeth.4402) and DPT [(Haghverdi *et al.*, Nat. Meth.,
  2016)](https://doi.org/10.1038/nmeth.3971)

<img src="./paul15/figures/draw_graph_fr.png" height="200"><img src="./paul15/figures/aga_graph.png" height="220">

#### Data from [Nestorowa, Hamey, *et al.* (2016)](http://doi.org/10.1182/blood-2016-05-716480)

In [*nestorowa16*](nestorowa16), we analyze data for early hematopoietic differentation.

<img src="./nestorowa16/figures/draw_graph_fr.png" height="200"><img src="./nestorowa16/figures/aga_graph.png" height="220">

### Lineage tree for whole cell atlas of an adult animal

In [planaria](planaria), we reconstruct the lineage tree of the whole cell atlas
of planaria (Plass, Jordi *et al.*, in preparation, 2017).

<img src="./planaria/figures/tsne.png" height="350"><img src="./planaria/figures/aga_graph.svg" height="450">

### PBMC cells

All of the following datasets contain correct motifs of PBMC differentiation, but as they are quite strongly clustering, the global tree cannot be inferred.

#### 3K cells from [10X (2017)](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k)

<img src="./figs_pbmc3k/tsne.png" height="200"><img src="./figs_pbmc3k/aga_graph.png" height="220">

For details, see [here](pbmc3k.html).

#### 33K cells from [10X (2017)](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc33k)

<img src="./figs_pbmc33k/tsne.png" height="200"><img src="./figs_pbmc33k/aga_graph.png" height="220">

For details, see [here](pbmc33k.html).

#### 68K cells from [Zheng *et al.* (2017)](https://doi.org/10.1038/ncomms14049)

<img src="./figs_zheng17/tsne.png" height="200"><img src="./figs_zheng17/aga_graph.png" height="220">

For details, see [here](zheng17.html).





