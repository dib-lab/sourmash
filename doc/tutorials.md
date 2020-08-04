# sourmash tutorials and notebooks

## The first two tutorials!

These tutorials are both command line tutorials that should work on Mac OS
X and Linux. They require about 5 GB of disk space and 5 GB of RAM.

* [The first sourmash tutorial - making signatures, comparing, and searching](tutorial-basic.md)

* [Using sourmash LCA to do taxonomic classification](tutorials-lca.md)

## Background and details

These next three tutorials are all notebooks that you can view, run
yourself, or run interactively online via the
[binder](http://mybinder.org) service.

* [An introduction to k-mers for genome comparison and analysis](kmers-and-minhash.md)

* [Some sourmash command line examples!](sourmash-examples.md)

* [Working with private collections of signatures.](sourmash-collections.md)

## More information

For more information on sequencing data with sourmash, check out our [ANGUS tutorial.](angus-tutorial.md)

If you are a Python programmer, you might also be interested in our [API examples](api-example.md) as well as a short guide to [Using the `LCA_Database` API.](using-LCA-database-API.ipynb)

If you prefer R, we have [a short guide to using sourmash output with R](other-languages.md).

## Customizing matrix and dendrogram plots in Python

If you're interested in customizing the output of `sourmash plot`,
which produces comparison matrices and dendrograms, please see
[Building plots from `sourmash compare` output](plotting-compare.md).

## Contents:

```{toctree}
:maxdepth: 2

tutorial-basic
tutorials-lca
kmers-and-minhash
sourmash-examples
sourmash-collections
angus-tutorial
api-example
using-LCA-database-API
other-languages
```
