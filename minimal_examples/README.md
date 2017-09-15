*First compiled: September 13, 2017.*

## Studying simple and clean simulated datasets

### Contents

* [*sim_data*](sim_data.ipynb) contains the main results
* [*robustness*](robustness.ipynb) studies the robustness of these results
* [*comparisons*](comparisons) discusses results obtained with other algorithms
  ([*comparisons_exports*](comparisons_exports.ipynb) exports data)

### Summary Tree Inference

Here, we study a simple, almost noise-free simulated dataset, which serves as a
minimal example with an unambiguous ground truth
([*krumsiek11_blobs/X.csv*](comparisons/data/krumsiek11_blobs/X.csv)). The data
contains two clusters and a continous tree-like manifold associated with
simulated hematopoietic differentiation. Below, we reproduce the abstracted
graph and the inferred tree of Figure 1 of the paper:

<img src="./figures/aga_draw_graph_fr.png" height="200"><img src="./figures/aga_graph.png" height="200">

### Comparison with competing algorithms

In [*comparisons*](comparisons), we show that no competing algorithm
is able to produce a result for the minimal dataset that is even close to
correct.

### Robustness

If we change parameters in graph abstraction, we obtain seemingly very different
abstractions of the data. Here, for a fine-grained resolution

<img src="./figures/aga.png" height="200">

and here, with a very coarse-grained resolution:

<img src="./figures/aga_new.png" height="200">

Upon closer inspection, one realizes that all results represent the same
topology and are correct. To measure how much two tree topologies differ, we
suggest to compare all paths between leaf nodes in the inferred abstracted
graphs. For example,
```
      path = ['21', '8', '18', '7', '9', '2'],
path_mapped = [['7', '2'], ['6', '7', '2'], ['2', '7'], ['2', '9'], ['9', '10', '3'], ['11', '10']],
      path_new = ['7', '2', '9', '10', '11'],
-> n_agreeing_steps = 4 / n_steps = 4.
```
agree.

For the simulated dataset, graph abstraction almost always yields the correct
topology; no matter which parameters are chosen.

<img src="./figures/summary.png" height="230">

See [*robustness*](robustness.ipynb) for all details.

### Different degrees of clustering

In [*sim_data*](sim_data.ipynb), we also show how graph abstraction behaves on
data with different degrees of clustering.

<img src="./figures/aga_cluster_std1.png" height="200">
<img src="./figures/aga_cluster_std6.png" height="200">
<img src="./figures/aga_cluster_std10.png" height="200">
