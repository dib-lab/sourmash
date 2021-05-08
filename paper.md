---
title: 'sourmash: a tool to quickly search, compare, and analyze genomic and metagenomic data sets'
tags:
  - MinHash
  - k-mers
  - Python
  - Rust
authors:
 - name: Luiz Irber
   orcid: 0000-0003-4371-9659
   affiliation: University of California, Davis
 - name: Tim Head
   orcid: 0000-0003-0931-3698
   affiliation: @@@
 - name: Harriet Alexander
   orcid: 0000-0003-1308-8008
   affiliation: Woods Hole Oceanographic Institute
 - name: Olga Botvinnik
   orcid: 0000-0003-4412-7970
   affiliation: Chan-Zuckerberg Biohub
 - name: Phillip Brooks
   orcid: 0000-0003-3987-244X
   affiliation: 10x Genomics
 - name: Laurent Gautier
   orcid: 0000-0003-0638-3391
   affiliation: @@@
 - name: Lisa Johnson
   orcid: 0000-0002-3600-7218
   affiliation: 10x Genomics
 - name: Fabian Klötzl
   orcid: @@@
   affiliation: @@@
 - name: David Koslicki
   orcid: 0000-0002-0640-954X
   affiliation: Pennsylvania State University
 - name: Katrin Leinweber
   orcid: 0000-0001-5135-5758
   affiliation: @@@
 - name: Ivan Ogasawara
   orcid: @@@
   affiliation: @@@
 - name: N. Tessa Pierce
   orcid: 0000-0002-2942-5331
   affiliation: University of California, Davis
 - name: Taylor Reiter
   orcid: 0000-0002-7388-421X
   affiliation: University of California, Davis
 - name: Camille Scott
   orcid: 0000-0001-8822-8779
   affiliation: University of California, Davis
 - name: Connor Skennerton
   orcid: 0000-0003-1320-4873
   affiliation: @@@
 - name: Daniel Standage
   orcid: 0000-0003-0342-8531.
   affiliation: @@@
 - name: Joshua Swamidass
   orcid: 0000-0003-2191-0778
   affiliation: @@@
 - name: Connor Tiffany
   orcid: @@@
   affiliation: @@@
 - name: Erik Young
   orcid: @@@
   affiliation: University of California, Davis
 - name: C. Titus Brown
   orcid: 0000-0001-6001-2677
   affiliation: University of California, Davis
date: 3 Mar 2021
bibliography: paper.bib
---

# Summary

sourmash is a command line tool and Python library for sketching
collections of DNA, RNA, and amino acid k-mers for biological sequence
search, comparison, and analysis [@Pierce:2019]. The most recent
release, sourmash v4, is built on top of Rust and provides an
experimental Rust interface.

MinHash sketches provide a lightweight way to store "signatures" of
large DNA or RNA sequence collections, and then compare or search them
using a Jaccard index [@Ondov:2015].  MinHash sketches can be used
to identify samples, find similar samples, identify data sets with
shared sequences, and build phylogenetic trees.

Since sourmash v1 was released in 2016 [@Brown:2016], sourmash has expanded
to support new database types and many more command line functions.
In particular, sourmash now has robust support for both Jaccard similarity
and containment calculations, which enables analysis and comparison of data sets
of different sizes, including large metagenomic samples.

# Statement of Need

Large collections of genomes, transcriptomes, and raw sequencing data
sets are readily available in biology, and the field needs lightweight
computational methods for searching and summarizing the content of
both public and private collections. sourmash provides a flexible set
of programmatic functionality for this purpose, together with a robust
and well-tested command-line interface. It has been used in over 70
publications (based on citations of @Brown:2016) and it continues
to expand in functionality.

# Acknowledgements

This work is funded in part by the Gordon and Betty Moore Foundation’s
Data-Driven Discovery Initiative [GBMF4551 to CTB].

# References
