## ReScience Ten Year Reproducibility Challenge

[![DOI](https://zenodo.org/badge/238372063.svg)](https://zenodo.org/badge/latestdoi/238372063)

This repository is my contribution to the
[ReScience Ten Years Reproducibility
Challenge](http://rescience.github.io/ten-years).

> Broman KW (2020) Reproducibility report: Identifying essential genes by mutagenesis. ReScience C 6 (1): #12
> [doi:10.5281/zenodo.3959516](https://doi.org/10.5281/zenodo.3959516)

I sought to reproduce the following paper.

> Lamichhane G, Zignol M, Blades NJ, Geiman DE, Dougherty A, Broman KW, Bishai WR (2003)  A post-genomic method
> for predicting essential genes at subsaturation levels of mutagenesis:
> application to _Mycobacterium tuberculosis_.  [Proc Natl Acad Sci USA](http://www.pnas.org/) 100:7213-7218
> [![PubMed](https://kbroman.org/pages/icons16/pubmed-icon.png)](https://www.ncbi.nlm.nih.gov/pubmed/12775759)
> [![pdf (260k)](https://kbroman.org/pages/icons16/pdf-icon.png)](https://www.pnas.org/content/pnas/100/12/7213.full.pdf)
> [![R/negenes](https://kbroman.org/pages/icons16/R-icon.png)](https://github.com/kbroman/negenes)
> [![doi](https://kbroman.org/pages/icons16/doi-icon.png)](https://doi.org/10.1073/pnas.1231432100)

A snapshot of my original project directory for my work related to
that paper is in
a [separate repository on
github](https://github.com/kbroman/Project_Lamichhane2003).

The primary software is the R/negenes package for R, available on
[github](https://github.com/kbroman/negenes) and
[CRAN](https://cran.r-project.org/package=negenes).

In this repository:

- [`negenes/`](negenes) contains the [R/negenes](https://github.com/kbroman/negenes)
  package, version 1.0-12 (2019-08-05).

- [`original/`](original) contains a subset of the [original project
  directory](https://github.com/kbroman/Projects_Lamichhane2003). Much
  of the material that was not used for the paper
  was deleted, as were intermediate results (`.RData` files).

- [`talk/`](talk) contains a subset of the [repository for slides for
  a 2002 talk on the work](https://github.com/kbroman/Talk_Mtb);
  deleted the figures, slides PDF, and a big .RData file. The key
  thing here is [`talk/R/circlefig.R`](talk/R/circlefig.R) which was
  used to create Figure 1b (of the circular genome).

- [`reproduction/`](reproduction) contains my reproduction of the
  analysis, as an [R Markdown
  document](reproduction/reproduction.Rmd).
  The compiled document can be [viewed on the
  web](https://kbroman.org/Paper_ReScience2020/reproduction/reproduction.html).

- [`article/`](article) contains the article about the reproduction effort.
  The [article PDF can be viewed on the web](https://kbroman.org/Paper_ReScience2020/article/article.pdf).


### Dependencies

The reproduction of the analysis, in [`reproduction`](reproduction),
relies on the following tools:

- [R](https://www.r-project.org)

- R packages: [negenes](https://github.com/kbroman/negenes),
  [rmarkdown](https://github.com/rstudio/rmarkdown),
  [devtools](https://devtools.r-lib.org/),
  [R.utils](https://github.com/HenrikBengtsson/R.utils),
  [xtable](http://xtable.r-forge.r-project.org/),
  [readxl](https://readxl.tidyverse.org/),
  [gt](https://gt.rstudio.com)

  Install the first four with `install.packages()`.

  ```r
  install.packages(c("negenes", "rmarkdown", "devtools", "R.utils", "xtable", "readxl", "broman"))
  ```

  Then install [gt](https://gt.rstudio.com) with `devtools::install_github()`.

  ```r
  devtools::install_github("rstudio/gt")
  ```

- [Perl](https://www.perl.org/)

- [Pandoc](https://pandoc.org/) (most easily installed by installing
  [RStudio Desktop](https://rstudio.com/products/rstudio/download/#download))

- [GNU Make](https://www.gnu.org/make)


### License

This is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, is available at https://www.r-project.org/Licenses/GPL-3
