[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![last commit](https://img.shields.io/github/last-commit/ocbe-uio/rBAPS.svg)](https://github.com/ocbe-uio/rBAPS/commits/master)
[![code size](https://img.shields.io/github/languages/code-size/ocbe-uio/rBAPS.svg)](https://github.com/ocbe-uio/rBAPS)
[![check status](https://github.com/ocbe-uio/rBAPS/workflows/R-CMD-check/badge.svg)](https://github.com/ocbe-uio/rBAPS/actions)
[![codecov](https://codecov.io/gh/ocbe-uio/rBAPS/branch/develop/graph/badge.svg)](https://codecov.io/gh/ocbe-uio/rBAPS)
[![CodeFactor](https://www.codefactor.io/repository/github/ocbe-uio/rBAPS/badge)](https://www.codefactor.io/repository/github/ocbe-uio/rBAPS)

# rBAPS
R implementation of the compiled Matlab BAPS software for Bayesian Analysis of Population Structure.

## Installation

rBAPS is currently under development and a stable version is yet to be released. You can install the development version by issuing the command below on an interactive R session.

```r
remotes::install_github("ocbe-uio/rBAPS", "develop")
```

At the current stage, we still recommend using the original BAPS software for production.

## Contributing

rBAPS is Open Software licenced under the [GPL-3](https://tldrlegal.com/license/gnu-general-public-license-v3-(gpl-3)), so all contributions are welcome. You can find a list of tasks we're working on for the first release of rBAPS [here](https://github.com/ocbe-uio/rBAPS/milestone/1). and pitfalls on the [PITFALLS.md](PITFALLS.md) file.

Found something not to be working as you think it should? Please help us by [opening an issue](https://github.com/ocbe-uio/rBAPS/issues/new).

## References

### How to cite this project

After installing the package, issue `citation("rBAPS")` in R. You will get a pure text and a BibTeX entry of rBAPS. Alternatively, you can use the [CITATION.cff](CITATION.cff) file, which conforms to the [Citation File Format](https://citation-file-format.github.io/) standard.

### Scientific papers

These are some of the scientific publications that served as base for the creation of rBAPS:

- Corander, J. and Marttinen, P. (2006), Bayesian identification of admixture events using multilocus molecular markers. Molecular Ecology, 15: 2833-2843. doi:[10.1111/j.1365-294X.2006.02994.x](https://doi.org/10.1111/j.1365-294X.2006.02994.x)
- Corander, J., Sirén, J., & Arjas, E. (2008). Bayesian spatial modeling of genetic population structure. Computational Statistics, 23(1), 111-129. doi:[10.1007/s00180-007-0072-x](https://link.springer.com/content/pdf/10.1007/s00180-007-0072-x.pdf)
- Corander, J., Marttinen, P., Sirén, J. et al. Enhanced Bayesian modelling in BAPS software for learning genetic structures of populations. BMC Bioinformatics 9, 539 (2008) doi:[10.1186/1471-2105-9-539](https://doi.org/10.1186/1471-2105-9-539)
- Tonkin-Hill G, Lees JA, Bentley SD et al. RhierBAPS: An R implementation of the population clustering algorithm hierBAPS [version 1; peer review: 2 approved]. Wellcome Open Res 2018, 3:93 (https://doi.org/10.12688/wellcomeopenres.14694.1)

### Software

Below we list some of the software directly related to rBAPS:

- [Original BAPS software](http://www.helsinki.fi/bsg/software/BAPS/)
- [RhierBAPS: An R implementation of the population clustering algorithm hierBAPS](https://github.com/gtonkinhill/rhierbaps)
- [fastbaps: A fast genetic clustering algorithm that approximates a Dirichlet Process Mixture model](https://github.com/gtonkinhill/fastbaps)
