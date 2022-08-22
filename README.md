# Contributions of adaptation and purifying selection to SARS-CoV-2 evolution

__Richard A. Neher__


Continued evolution and adaptation of SARS-CoV-2 has lead to more transmissible and immune-evasive variants with profound impact on the course of the pandemic.
Here I analyze the evolution of the virus over 2.5 years since its emergence and estimate rates of evolution for synonymous and non-synonymous changes separately for evolution within clades -- well defined mono-phyletic groups with gradual evolution -- and for the pandemic overall.
The rate of synonymous mutations is found to be around 6 changes per year.
Synonymous rates within variants vary little from variant to variant and are compatible with the overall rate.
In contrast, the rate at which variants accumulate amino acid changes (non-synonymous mutation) was initially around 12-16 changes per year, but in 2021 and 2022 dropped to 6-9 changes per year.
The overall rate of non-synonymous evolution, that is across variants, is estimated to be about 25 amino acid changes per year.
This 2-fold higher rate indicates that the evolutionary process that gave rise to the different variants is qualitatively different from that in typical transmission chains and likely dominated by adaptive evolution.
I further quantify the spectrum of mutations and purifying selection in different SARS-CoV-2 proteins.
Many accessory proteins evolve under limited evolutionary constraint with little short term purifying selection.
About half of the mutations in other proteins are strongly deleterious and rarely observed, not even at low frequency.


### Repository structure

This repository contains scripts and source files associated with a manuscript on SARS-CoV-2 virus evolution.

The analysis can be run using `snakemake` and requires standard python libraries, as well as `treetime` (`phylo-treetime` in pip).

The analysis can be run using the open data files provisioned by Nextstrain (download default in the workflow).

The directory `manuscript` contains the `TeX` files associated with manuscript, the bibliography, and the figures.

The `data` directory contains some derived files, including the rate estimates, the mutation distribution, and fitness cost landscapes.

