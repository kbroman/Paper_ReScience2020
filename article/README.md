### Article about reproducing the analysis in Lamichhane et al. (2003)

This repository contains the article about the reproduction work,
using the [ReScience LaTeX template](https://github.com/rescience/template).

The key customized files are:

- [`metadata.yaml`](metadata.yaml)
- [`bibliography.bib`](bibliography.bib)
- [`content.tex`](content.tex)

The figures and tables are drawn in from directories above: the
original figures in [`../original/`](original/) and [`../talk/`](talk/),
the reproduced figures in
[`../reproduction/Figs/`](../reproduction/Figs/), and the LaTeX tables
in [`../reproduction/Tabs/`](../reproduction/Tabs/).

The compiled document is in [`article.pdf`](article.pdf),
also available at
<https://kbroman.org/Paper_ReScience2020/article/article.pdf>.

To compile the article type `make`.
