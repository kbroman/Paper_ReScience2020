## ReScience Ten Year Reproducibility Challenge

This repository is my contribution to the
[ReScience Ten Years Reproducibility
Challenge](http://rescience.github.io/ten-years). I'm seeking to
reproduce the following paper.

> Lamichhane G, Zignol M, Blades NJ, Geiman DE, Dougherty A, **Broman KW**, Bishai WR (2003)  A post-genomic method
> for predicting essential genes at subsaturation levels of mutagenesis:
> application to <i>Mycobacterium tuberculosis</i>.  [Proc Natl Acad Sci USA](http://www.pnas.org/) 100:7213-7218
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
