@author Till Budde (tilljanis.budde@tu-dortmund.de)
@date 16.09.2021

## Introduction

Within the field of computational mechanics the implementation of material
interfaces has been of interest for a while.
At it's most basic material interfaces exist at any point where the surrounding
material exhibits some form of discontinuity, wether this is the discontinuity of
displacements as seen in cracks and shear bands, or the discontinuity of material
as is the case with the surface of enclosures depends of the specific problem
statement. 
A generalized framework for material interfaces can therefore proof useful
in a wide area of study, from polycristalyne materials to the calculation of 
adhesive bonds.
In this work a new finite-element Code, suited for the usage of interface-elements,
is build on the basis of the finite-element library dealii.
While dealii is a powerful and flexible framework for fem calculations, it is 
a far cry from a commercial user-centric fem software.
As such a considerable amount of work was spent in writing the basic fem code,
with it's many adaptations specific to the use of interface-elements.
Furthermore dealii has a quite stringent design philosophy which often clashes
with the needs of a generalized interface framework, starting with the creation
of the mesh (and its zero volume interface elements), to the usage of the 
"masterelement" class FEFaceValues and the assembly of the system of equations.
Nevertheless the creation of the FEM code was mostly succesful.
On the basis of this framework a library of Interface elements with differing 
material models could be implemented.

## Goals of this work

- [ ] Implementation of a FEM Code suited to the usage of interface elements
- [ ] Implementation of a library of interface elements with different material models
  - [X] Spring like interface
  - [X] Damaged Fiber Model interface 
  - [ ] Elastoplastic fiber model
