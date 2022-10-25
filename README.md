# neuronalActionPotential

### About

neuronalActionPotential is a C++ code for modeling the action potential propagation along the neuronal axon. The branch cableHH models the cable theory based one dimensional Hodgkin-Huxley electrical network model. The classical work of Hodgkin and Huxley was a landmark model in terms of providing deep insights into the ionic basis of action potential propagation in nerve cells. The primal field that is solved for is the voltage. This model is applied to numerically estimate the conduction velocity in a rat axon. This electrical network model takes into account the membrane capacitance and the
ionic currents due to the sodium ions, potassium ions and some leak current through the respective ion channels in the neuronal membrane.

### Installation:

neuronalActionPotential code builds on top of the deal.II library.

1) Install CMake, PETSc, Trilinos, SLEPc, p4est, and deal.II (version 9.3.0 recommended)<br>

2) Clone the neuronalActionPotential GitHub repository <br>
```
$ git clone https://github.com/cmmg/neuronalActionPotential.git
$ cd neuronalActionPotential
$ git checkout cableHH
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
