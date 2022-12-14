# neuronalActionPotential

### About

neuronalActionPotential is a C++ code for modeling the action potential propagation along the neuronal axon. The branch pnpPeriaxonal models the electro-diffusive PNP considering the presence of the myelin sheath and the submyelin peri-axonal space in a rat axon. 

This framework for modeling spatio-temporal ionic charge distribution and the resulting voltage evolution using high-fidelity partial differential equations (PDE) modeling coupled electrostatics and electrochemistry can be capable of faithfully representing spatio-temporally heterogenous evolution of the ionic and voltage distributions leading to generation, propagation and potentially disruption of the neuronal action potential. The primal fields that are solved for are the voltage and the concentration of sodium/ potassium/ chloride ions. This electro-diffusive model is applied to numerically estimate the conduction velocity in a rat axon. The salient features of the computational implementation are: adaptive mesh refinement near the nodes of Ranvier, adaptive time-stepping schemes, support for parallel direct and iterative (Krylov-subspace) solvers with Jacobi/SOR preconditioning.


### Installation:

neuronalActionPotential code builds on top of the deal.II library.

1) Install CMake, PETSc, Trilinos, SLEPc, p4est, and deal.II (version 9.3.0 recommended)<br>

2) Clone the neuronalActionPotential GitHub repository <br>
```
$ git clone https://github.com/cmmg/neuronalActionPotential.git
$ cd neuronalActionPotential
$ git checkout pnpPeriaxonal
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
