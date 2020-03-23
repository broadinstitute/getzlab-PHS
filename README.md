# Passenger Hotspots
Code to reproduce the figures and results from the manuscript "Passenger Hotspot Mutations in Cancer" (Hess _et al._ 2019, https://doi.org/10.1016/j.ccell.2019.08.002)

### What's in this repo?

This repository contains:
* MATLAB notebooks to reproduce all of the analyses from the manuscript in [`PHS/`](https://github.com/broadinstitute/getzlab-PHS/tree/master/PHS)
* Functions specific to this manuscript invoked by notebooks in [`misc_functions/`](https://github.com/broadinstitute/getzlab-PHS/tree/master/misc_functions)
* Code required to compile and run the [log-normal-Poisson regression model](https://github.com/broadinstitute/getzlab-LNP/) in [`LNP/`](https://github.com/broadinstitute/getzlab-LNP/tree/0bf0abe540f37d6b243fad05ad68dc7dee5a8deb)
* Code required for these analyses but not specific to this manuscript in [`funcs/`](https://github.com/broadinstitute/getzlab-LNPHS_utils/tree/296448f9de9a3079d39265470915c3c658162cd8)

### What's *not* in this repo?

* Any reference data used by these analyses, which totals approximately 19 GB.  These data are hosted in a Google storage bucket (`gs://getzlab-passengerhotspots`); to download this, please [install gsutil](https://cloud.google.com/storage/docs/gsutil_install) and run (in the root directory of this repo):

```
gsutil -m cp -r gs://getzlab-passengerhotspots/* ref
```

* Any dbGaP protected data, which will have to be obtained by individuals with appropriate authorization.

All external reference data are assumed by the code to reside in the `ref/` directory, and external mutation data in the `mutation_data/` directory. Protected data, when present, is clearly denoted in the code.

### How to run this code

All notebooks are intended to be run interactively in the MATLAB console.  Start MATLAB (any version R2014b or newer should work) in the root directory of this repo — this is necessary for `startup.m` to properly add dependencies to the MATLAB path. Code tested to work only under 64 bit Linux; other architectures may work after recompiling C/C++ `.mex` files.

## What's in each notebook?

The notebooks are roughly sorted in the order in which they should be run. However, grouping similar code together takes precedence over having a perfectly linear notebook structure.

Sections of the code that produce figures and tables are clearly denoted in the code as, e.g., "`Figure 3A`", "`Figure S5C`", or "`Table S1`".

* [`00_process_maf`](https://github.com/broadinstitute/getzlab-PHS/blob/master/PHS/00_process_maf.m)
  * Take the raw mutations calls (in TSV format) and convert them to a format optimized for downstream analysis.
  * Filter regions of low mappability; apply addition panel-of-normals filtering
  * Split into subcohorts whose patients' tumors are dominated by specific mutational processes
  * **Generate Figure S4**, plots of the trinucleotide context distributions for mutations in these subcohorts

* [`01_process_territory`](https://github.com/broadinstitute/getzlab-PHS/blob/master/PHS/01_process_territory.m)
  * Calculate regions of low mappability
  * Process regression covariates
  * Tabulate sequence context territories

* [`02_run_sig_algos`](https://github.com/broadinstitute/getzlab-PHS/blob/master/PHS/02_run_sig_algos.m)
  * Run each of the four significance methods on the mutation dataset
* [`05_compute_LNP_posterior_predictives`](https://github.com/broadinstitute/getzlab-PHS/blob/master/PHS/05_compute_LNP_posterior_predictives.m)
  * Compute _p_-values for LNP regression
  * **Generate Figure S5C**, contrasting the effect of incorporating APOBEC3A substrate optimality ([Buisson _et al._ 2019](https://science.sciencemag.org/content/364/6447/eaaw2872)) on LNP _p_-values
* [`10_sig_analysis_and_q_scatter`](https://github.com/broadinstitute/getzlab-PHS/blob/master/PHS/10_sig_analysis_and_q_scatter.m)
  * **Generate Figure S2A**, scatterplots contrasting methods' _p_-values
* [`12_histogram_figure`](https://github.com/broadinstitute/getzlab-PHS/blob/master/PHS/12_histogram_figure.m)
  * **Generate Figure 3A**, the observed fraction of synonymous hotspots and expected fractions predicted by each model
  * **Generate Figure 3B**, visualizing probabilities of specific driver/passenger mutations as predicted by each model
* [`20_tabulate_effect_territories`](https://github.com/broadinstitute/getzlab-PHS/blob/master/PHS/20_tabulate_effect_territories.m)
  * Tabulate joint sequence context/protein coding effect territories (e.g., number of T(C->G)A mutations that cause synonymous mutations)
* [`21_global_effect_analysis`](https://github.com/broadinstitute/getzlab-PHS/blob/master/PHS/21_global_effect_analysis.m)
  * Compute expected fraction of protein coding effects given a set of mutations, if those mutations were randomly distributed throughout the exome
  * **Generate Figure 1**, contrasting the observed/expected fraction of protein coding effects for hotspot mutations significant by each model
  * **Generate Figure S1**, showing the distribution of protein coding effects for each trinucleotide context
* [`22_gene_dNdS_analysis`](https://github.com/broadinstitute/getzlab-PHS/blob/master/PHS/22_gene_dNdS_analysis.m)
  * Compute somatic dN/dS for each gene (for true/false positive truth sets)
  * **Generate Table S1**
  * **Generate Table S3**
* [`23_ROC-FDR`](https://github.com/broadinstitute/getzlab-PHS/blob/master/PHS/23_ROC-FDR.m)
  * **Generate Figure 2A**, receiver-operator characteristic analysis of the four methods
  * **Generate Figure 2B**, empirical FDR analysis of the methods
  * **Generate Figure S2B**, ROC plots using an alternate true positive truth set
* [`24_qq`](https://github.com/broadinstitute/getzlab-PHS/blob/master/PHS/24_qq.m)
  * **Generate Figure S2C**, quantile-quantile plots of the four methods' _p_-values
* [`30_hypermut`](https://github.com/broadinstitute/getzlab-PHS/blob/master/PHS/30_hypermut.m)
  * **Generate Figure S5A**, contrasting properties of mutational processes in hypermutated vs. non-hypermutated samples
* [`50_signature_analysis`](https://github.com/broadinstitute/getzlab-PHS/blob/master/PHS/50_signature_analysis.m)
  * **Generate Figure 4A**, visualizing log-normal-Poisson posterior distributions for different mutations processes
  * **Generate Figure 4B**, computing variance explained by each genomic covariate for different mutations processes
  * **Generate Figure S5B**, computing the variance explained by adding XR-seq coverage as a covariate to the UV mutational process
