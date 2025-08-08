spareGParts: Miscellaneous Emulation Methods
================

[![License: BSD
3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![](https://img.shields.io/badge/devel%20version-0.1.0-purple.svg)](https://github.com/knrumsey/khaos)

<!-- README.md is generated from README.Rmd. Please edit that file -->

### Installation

To install this package, use

``` r
# install.packages("devtools")
devtools::install_github("knrumsey/spareGParts")
```

The `README.Rmd` file can be knit locally for proper rendering of the
equations. Will add a link here if I ever write a short paper for this.

### Description

The `spareGParts` R package implements and/or packages several emulation
techniques (mostly approximate GPs) that are discussed in the literature
but without a convenient R implementation. Currently includes:

1.  `mpgp`: **Matching Pursuit GP**. A subset of data approach proposed
    by Keerthi and Chu (2005) with complexity O(nm^2).
2.  `svechgp` **Scaled Vecchia GP**. The scaled Vecchia GP approach of
    Katzfuss et al. (2020). *Note: I am not the author of this function.
    See [here](https://github.com/katzfuss-group/scaledVecchia) for
    original implementation. I am simply packagizing the code for
    convenience with permission from Matthias Katzfuss*.
3.  `rvm` **Relevance Vector Machine**. An implementation of the RVM
    algorithm via Tipping 2001. The `kernlab` implementation is faster
    (written in C++) but fails to yield probabilistic modeling. This
    implementation allows allows for specifying a discrete prior over
    the lengthscale parameter (and admits parallelization via
    `mclapply`). A pre-screening LASSO step is also available to bound
    the complexity of the algorithm to O(maxiter \* max_basis^3).
    Posterior probabilities are given by
    $$p(\ell_j | {\bf X, y}) \propto p(\ell_j) \max_{{\bf \alpha}, \sigma^2}p({\bf y} | \ell_j, {\bf\alpha}, \sigma^2)$$
    where $p({\bf y} | \ell_j, {\bf\alpha}, \sigma^2)$ is given by
    equation 7 in Tipping (2001).

### References

Tipping, Michael E. “Sparse Bayesian learning and the relevance vector
machine.” Journal of machine learning research 1.Jun (2001): 211-244.

Keerthi, Sathiya, and Wei Chu. “A matching pursuit approach to sparse
gaussian process regression.” Advances in neural information processing
systems 18 (2005).

Liu, Haitao, et al. “When Gaussian process meets big data: A review of
scalable GPs.” IEEE transactions on neural networks and learning systems
31.11 (2020): 4405-4423.

Chalupka, Krzysztof, Christopher KI Williams, and Iain Murray. “A
framework for evaluating approximation methods for Gaussian process
regression.” The Journal of Machine Learning Research 14.1 (2013):
333-350.

Ranjan, P., Haynes, R., and Karsten, R. (2011). A Computationally Stable
Approach to Gaussian Process Interpolation of Deterministic Computer
Simulation Data, Technometrics, 53(4), 366 - 378.

Katzfuss, Matthias, Joseph Guinness, and Earl Lawrence. “Scaled Vecchia
approximation for fast computer-model emulation.” SIAM/ASA Journal on
Uncertainty Quantification 10.2 (2022): 537-554.

### Copyright

Reference ID: \#O4792

*© 2024. Triad National Security, LLC. All rights reserved. This program
was produced under U.S. Government contract 89233218CNA000001 for Los
Alamos National Laboratory (LANL), which is operated by Triad National
Security, LLC for the U.S. Department of Energy/National Nuclear
Security Administration. All rights in the program are reserved by Triad
National Security, LLC, and the U.S. Department of Energy/National
Nuclear Security Administration. The Government is granted for itself
and others acting on its behalf a nonexclusive, paid-up, irrevocable
worldwide license in this material to reproduce, prepare. derivative
works, distribute copies to the public, perform publicly and display
publicly, and to permit others to do so.*
