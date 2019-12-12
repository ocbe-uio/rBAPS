# rBAPS
R implementation of the compiled Matlab BAPS software for Bayesian Analysis of Population Structure.

## Installation

rBAPS is currently under development and a stable version is yet to be released. You can install the development version by issuing the command below on an interactive R session.

```{r}
remotes::install_github("ocbe-uio/rBAPS", "dev")
```

## References

### Scientific papers

- Corander, J. and Marttinen, P. (2006), Bayesian identification of admixture events using multilocus molecular markers. Molecular Ecology, 15: 2833-2843. doi:[10.1111/j.1365-294X.2006.02994.x](https://doi.org/10.1111/j.1365-294X.2006.02994.x)
- Corander, J., Sirén, J., & Arjas, E. (2008). Bayesian spatial modeling of genetic population structure. Computational Statistics, 23(1), 111-129. doi:[10.1007/s00180-007-0072-x](https://link.springer.com/content/pdf/10.1007/s00180-007-0072-x.pdf)
- Corander, J., Marttinen, P., Sirén, J. et al. Enhanced Bayesian modelling in BAPS software for learning genetic structures of populations. BMC Bioinformatics 9, 539 (2008) doi:[10.1186/1471-2105-9-539](https://doi.org/10.1186/1471-2105-9-539)
- Tonkin-Hill G, Lees JA, Bentley SD et al. RhierBAPS: An R implementation of the population clustering algorithm hierBAPS [version 1; peer review: 2 approved]. Wellcome Open Res 2018, 3:93 (https://doi.org/10.12688/wellcomeopenres.14694.1)

### Software

- [Original BAPS software](http://www.helsinki.fi/bsg/software/BAPS/)
- [RhierBAPS: An R implementation of the population clustering algorithm hierBAPS](https://github.com/gtonkinhill/rhierbaps)
- [fastbaps: A fast genetic clustering algorithm that approximates a Dirichlet Process Mixture model](https://github.com/gtonkinhill/fastbaps)