# FLARE - Finite eLement Analysis foR Electromagnetics
## Quick Start
The `run_FLARE.m` file contains a few examples involving the different solvers in FLARE.
### Generating mesh files
Meshing in FLARE is performed using the open source mesh generator `gmsh` ([download gmsh here](https://gmsh.info/#Download)).

Once `gmsh` in installed in your system, run one of the .geo files in the `mesh` directory. A `.m` file will be generated, that can be loaded by FLARE and used as a mesh for the given case.
### Source terms
The `sources` directory contains a few examples of source terms for physical problems solved by FLARE.
## Feedback
Feedback on FLARE may be sent to <arturo.popoli@unibo.it>
## Credits
FLARE is developed and mantained by Arturo Popoli, Andrea Cristofolini, Giacomo Pierotti and Fabio Ragazzi, at the University of Bologna - IT
## Aknowledgements
The inverse-Laplace solver in FLARE uses a Matlab implementation the following routine, by prof. L. D'Amore and colleagues:
> Luisa D'Amore, Guiliano Laccetti, and Almerico Murli. 1999. An implementation of a Fourier series method for the numerical inversion of the Laplace transform. ACM Trans. Math. Softw. 25, 3 (Sept. 1999), 279–305. https://doi.org/10.1145/326147.326148
```bibtex
@article{10.1145/326147.326148,
    author = {D'Amore, Luisa and Laccetti, Guiliano and Murli, Almerico},
    title = {An Implementation of a Fourier Series Method for the Numerical Inversion of the Laplace Transform},
    year = {1999},
    issue_date = {Sept. 1999},
    publisher = {Association for Computing Machinery},
    address = {New York, NY, USA},
    volume = {25},
    number = {3},
    issn = {0098-3500},
    url = {https://doi.org/10.1145/326147.326148},
    doi = {10.1145/326147.326148},
    journal = {ACM Trans. Math. Softw.},
    month = {sep},
    pages = {279–305},
    numpages = {27},
}
```
