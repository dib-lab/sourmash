# Using `sourmash` output with R and other languages

Most of the sourmash shell commands output CSV files upon request;
these files have headers and are straightforward to load into R.
Below are some code snippets and links that might be useful.

## R code for working with compare output

(by Taylor Reiter)

`sourmash compare` can output matrices in a CSV format, which can
easily be read into R.  For example, if you download the Eschericia
signature collection as in
[the sourmash tutorial](tutorial-basic.md#make-and-search-a-database-quickly),
then the shell command

```shell
sourmash compare ecoli_many_sigs/*.sig --csv ecoli.cmp.csv
```

will output a file `ecoli.cmp.csv` that can be loaded into R like so:

```r

sourmash_comp_matrix <- read.csv("ecoli.cmp.csv")

# Label the rows
rownames(sourmash_comp_matrix) <- colnames(sourmash_comp_matrix)

# Transform for plotting
sourmash_comp_matrix <- as.matrix(sourmash_comp_matrix)

```

[See the output of plotting and clustering this matrix](doc/_static/ecoli-cmp.html),
produced by [this RMarkdown file](doc/_static/ecoli-cmp.Rmd).

You can download the `ecoli.cmp.csv` file itself [here](doc/_static/ecoli.cmp.csv).

