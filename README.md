# Computing Radially-Symmetric Solutions of the Ultra-Relativistic Euler Equations with Entropy-Stable Discontinuous Galerkin Methods

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/TODO.svg)](https://doi.org/TODO)

This repository contains information and code to reproduce the results presented in the
article
```bibtex
@online{thein2025computing,
  title={{C}omputing Radially-Symmetric Solutions of the
         Ultra-Relativistic {E}uler Equations with Entropy-Stable
         Discontinuous {G}alerkin Methods},
  author={Thein, Ferdinand and Ranocha, Hendrik},
  year={2025},
  month={TODO},
  eprint={TODO},
  eprinttype={arxiv},
  eprintclass={math.NA}
}
```

If you find these results useful, please cite the article mentioned above. If you
use the implementations provided here, please **also** cite this repository as
```bibtex
@misc{thein2025computingRepro,
  title={Reproducibility repository for
         "{C}omputing Radially-Symmetric Solutions of the
         Ultra-Relativistic {E}uler Equations with Entropy-Stable
         Discontinuous {G}alerkin Methods"},
  author={Thein, Ferdinand and Ranocha, Hendrik},
  year={2025},
  howpublished={\url{https://github.com/ranocha/2025_ultrarelativistic_euler}},
  doi={TODO}
}
```

## Abstract

The ultra-relativistic Euler equations describe gases in the relativistic case when the thermal energy dominates.
These equations for an ideal gas are given in terms of the pressure, the spatial part of the dimensionless four-velocity, and the particle density.
Kunik et al. (2024, https://doi.org/10.1016/j.jcp.2024.113330) proposed genuine multi-dimensional benchmark problems for the ultra-relativistic Euler equations.
In particular, they compared full two-dimensional discontinuous Galerkin simulations for radially symmetric problems with solutions computed using a specific one-dimensional scheme.
Of particular interest in the solutions are the formation of shock waves and a pressure blow-up.
In the present work we derive an entropy-stable flux for the ultra-relativistic Euler equations.
Therefore, we derive the main field (or entropy variables) and the corresponding potentials.
We then present the entropy-stable flux and conclude with simulation results for different test cases both in 2D and in 3D.




## Numerical experiments

To reproduce the numerical experiments presented in this article, you need
to install [Julia](https://julialang.org/). The numerical experiments presented
in this article were performed using Julia v1.10.9.

First, you need to download this repository, e.g., by cloning it with `git`
or by downloading an archive via the GitHub interface. Then, you need to start
Julia in the `code` directory of this repository and follow the instructions
described in the `README.md` file therein.


## Authors

- [Ferdinand Thein](https://www.ferdinandthein.de) (Johannes Gutenberg University Mainz, Germany)
- [Hendrik Ranocha](https://ranocha.de) (Johannes Gutenberg University Mainz, Germany)


## License

The code in this repository is published under the MIT license, see the
`LICENSE` file.


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
