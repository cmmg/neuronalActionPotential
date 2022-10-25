# neuronalActionPotential

### About

neuronalActionPotential is a C++ code for modeling the action potential propagation along the neuronal axon. The branch cableDC models the double-cable electrical network model. The presence of the submyelin peri-axonal region has been proposed as a potential pathway for rapid electrical
conduction along the axon and is considered in the double-cable electrical circuit.

### Installation:

neuronalActionPotential code builds on top of the deal.II library.

1) Install CMake, PETSc, Trilinos, SLEPc, p4est, and deal.II (version 9.3.0 recommended)<br>

2) Clone the neuronalActionPotential GitHub repository <br>
```
$ git clone https://github.com/cmmg/neuronalActionPotential.git
$ cd neuronalActionPotential
$ git checkout cableDC
$ cmake .
$ make -j nprocs
  ```
[here nprocs denotes the number of processors]

### Visualization:

  Output of the primary fields is in the vtk format (parallel:*.pvtu, serial:*.vtu files). These can be visualized with the following open source applications:
  1. VisIt (https://visit-dav.github.io/visit-website/releases-as-tables/)
  2. Paraview (http://www.paraview.org/download/)


License
-------
GNU Lesser General Public License (LGPL). Please see the file LICENSE for details.
