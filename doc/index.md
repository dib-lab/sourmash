# Welcome to sourmash!

sourmash is a command-line tool and Python library for computing
[hash sketches](https://en.wikipedia.org/wiki/MinHash) from DNA
sequences, comparing them to each other, and plotting the results.
This allows you to estimate sequence similarity between even very
large data sets quickly and accurately.

sourmash can be used to quickly search large databases of genomes
for matches to query genomes and metagenomes; see [our list of
available databases](databases.md).

sourmash also includes k-mer based taxonomic exploration and
classification routines for genome and metagenome analysis. These
routines can use the NCBI and GTDB taxonomies but do not depend on them
specifically.

We have [several tutorials](tutorials.md) available! Start with
[Making signatures, comparing, and searching](tutorial-basic.md).

The sourmash software can be cited as
[Brown & Irber, 2019](https://doi.org/10.12688/f1000research.19675.1). The
paper
[Large-scale sequence comparisons with sourmash (Pierce et al., 2019)](https://doi.org/10.12688/f1000research.19675.1)
gives an overview of how sourmash works and what its major use cases
are.  Please also see the `mash`
[software](http://mash.readthedocs.io/en/latest/) and
[paper (Ondov et al., 2016)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x)
for background information on how and why MinHash works.

**Questions? Thoughts?** Ask us on the [sourmash issue tracker](https://github.com/dib-lab/sourmash/issues/)!

**Want to migrate to sourmash v4?** sourmash v4 is now available, and
has a number of incompatibilites with v2 and v3. Please see 
[our migration guide](support.md#migrating-from-sourmash-v3-x-to-sourmash-v4-x)!

----

To use sourmash, you must be comfortable with the UNIX command line;
programmers may find the [Python library and API](api.md) useful as well.

If you use sourmash, please cite us! In the body of the paper, we request
that you cite the F1000Research paper on use cases:

>"We used sourmash to do x, y and z. (Pierce et al., 2019).

and in the Methods, we suggest citing the specfic version of the software,
together with a reference to the JOSS software paper:

>"We used sourmash vA.B.C (Brown and Irber, 2016).

The full citations are:

>Pierce NT, Irber L, Reiter T et al. Large-scale sequence comparisons
>with sourmash [version 1; peer review: 2 approved]. F1000Research
>2019, 8:1006 (https://doi.org/10.12688/f1000research.19675.1)

and

> Brown and Irber (2016),
> **sourmash: a library for MinHash sketching of DNA**.
> Journal of Open Source Software, 1(5), 27, [doi:10.21105/joss.00027](https://joss.theoj.org/papers/3d793c6e7db683bee7c03377a4a7f3c9)

## sourmash in brief

sourmash uses MinHash-style sketching to create "signatures", compressed
representations of DNA/RNA sequence.  These signatures can then
be stored, searched, explored, and taxonomically annotated.

* `sourmash` provides command line utilities for creating, comparing,
  and searching signatures, as well as plotting and clustering
  signatures by similarity (see [the command-line docs](command-line.md)).

* `sourmash` can **search very large collections of signatures** to find matches
  to a query.

* `sourmash` can also **identify parts of metagenomes that match known genomes**,
  and can **taxonomically classify genomes and metagenomes** against databases
  of known species.

* `sourmash` can be used to **search databases of public sequences**
  (e.g. all of GenBank) and can also be used to create and search databases
  of **private sequencing data**.

* `sourmash` supports saving, loading, and communication of signatures
  via [JSON](http://www.json.org/), a ~human-readable and editable format.

* `sourmash` also has a simple Python API for interacting with signatures,
  including support for online updating and querying of signatures
  (see [the API docs](api.md)).

* `sourmash` relies on an underlying Rust core for performance.

* `sourmash` is developed [on GitHub](https://github.com/dib-lab/sourmash)
  and is **freely and openly available** under the BSD 3-clause license.
  Please see [the README](https://github.com/dib-lab/sourmash/blob/latest/README.md)
  for more information on development, support, and contributing.

You can take a look at sourmash analyses on real data
[in a saved Jupyter notebook](https://github.com/dib-lab/sourmash/blob/latest/doc/sourmash-examples.ipynb),
and experiment with it yourself
[interactively in a Jupyter Notebook](https://mybinder.org/v2/gh/dib-lab/sourmash/latest?filepath=doc%2Fsourmash-examples.ipynb)
at [mybinder.org](http://mybinder.org).

## Installing sourmash

You can use pip:
```bash
$ pip install sourmash
```

or conda:
```bash
$ conda install -c conda-forge -c bioconda sourmash
```

Please see [the README file in github.com/dib-lab/sourmash](https://github.com/dib-lab/sourmash/blob/latest/README.md)
for more information.

## Memory and speed

sourmash has relatively small disk and memory requirements compared to
many other software programs used for genome search and taxonomic
classification.

`sourmash search` and `sourmash gather` can be used to search 100k
genbank microbial genomes ([using our prepared databases](databases.md)
with about 20 GB of disk and in under 1 GB of RAM.
Typically a search for a single genome takes about 30 seconds on a laptop.

`sourmash lca` can be used to search/classify against all genbank
microbial genomes with about 200 MB of disk space and about 10 GB of
RAM. Typically a metagenome classification takes about 1 minute on a
laptop.

## sourmash versioning

We support the use of sourmash in pipelines and applications
by communicating clearly about bug fixes, feature additions, and feature
changes. We use version numbers as follows:

* Major releases, like v4.0.0, may break backwards compatibility at
  the command line as well as top-level Python/Rust APIs.
* Minor releases, like v4.1.0, will remain backwards compatible but
  may introduce significant new features.
* Patch releases, like v4.1.1, are for minor bug fixes; full backwards
  compatibility is retained.

If you are relying on sourmash in a pipeline or application, we
suggest specifying your version requirements at the major release,
e.g. in conda you would specify `sourmash>=3,<4`.

See [the Versioning docs](support.md) for more information on what our
versioning policy means in detail, and how to migrate between major
versions!

## Limitations

**sourmash cannot find matches across large evolutionary distances.**

sourmash seems to work well to search and compare data sets for
nucleotide matches at the species and genus level, but does not have much
sensitivity beyond that.  (It seems to be particularly good at
strain-level analysis.)  You should use protein-based analyses
to do searches across larger evolutionary distances.

**sourmash signatures can be very large.**

We use a modification of the MinHash sketch approach that allows us
to search the contents of metagenomes and large genomes with no loss
of sensitivity, but there is a tradeoff: there is no guaranteed limit
to signature size when using 'scaled' signatures.

## Logo

The sourmash logo was designed by Stéfanie Fares Sabbag,
with feedback from Clara Barcelos,
Taylor Reiter and Luiz Irber.

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img
alt="Creative Commons License" style="border-width:0"
src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />

The logo
is licensed under a <a rel="license"
href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons
Attribution-ShareAlike 4.0 International License</a>.

## Contents:

```{toctree}
:maxdepth: 2

command-line
tutorials
using-sourmash-a-guide
classifying-signatures
databases
api
more-info
support
developer
```

# Indices and tables

* {ref}`genindex`
* {ref}`modindex`
* {ref}`search`
