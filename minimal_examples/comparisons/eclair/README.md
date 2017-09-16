# Comparison with ECLAIR [(Giecold *et al.*, Nucl. Acids Research, 2016)](https://doi.org/10.1093/nar/gkw452)

## Contents

* [*eclair_plots*](eclair_plots.ipynb) visualizes the results of running eclair
* [*logfile_run_X_krumsiek11_scaled.txt*](logfile_run_X_krumsiek11_scaled.txt) shows the logging output and the parameters that allow to reproduce the results discussed [here](../)
* [*X_krumsiek11_scaled.txt*](X_krumsiek11_scaled.txt) contains the scaled data matrix for the simple tree, [*X_krumsiek11.txt*](X_krumsiek11.txt) contains the equivalent unscaled matrix and [*X_krumsiek11_blobs.txt*](X_krumsiek11_blobs.txt) contains the data describing tree and clusters; all in tab-separated format as required by ECLAIR. These files have been generated in [*../../comparison_exports*](../../comparison_exports.ipynb), which also provides visualizations of the data.

## Notes

* See the discussion of the results [here](../).

* We acknowledge help by G. Giecold and S. P. Garcia, in parts by email, in parts here https://github.com/GGiecold/ECLAIR/issues/3, who adviced to scale the data matrix and provided default parameters.

* The algorithm could not handle the unnormalized data matrices *X_krumsiek11_blobs.txt* and *X_krumsiek11.txt*.
