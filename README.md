# phd-code

Code for Conway's Ph.D., for Mathematica Version&nbsp;12.

Everything should also work in Version&nbsp;11, except the following:

* Nonlinear PDEs / boundary conditions passed to `NDSolve` for
  `` NDSolve`FEM` ``.
  
* `Plot3D` where the domain is given in the form `Element[{x,y}, mesh]]`,
  where `mesh` is an `ElementMesh`.
  
  In Version&nbsp;11.0 this results in the error message
  > Plot3D::idomdim: {x,y} ∈ mesh does not have a valid dimension
    as a plotting domain.
  
  A workaround is to replace `mesh` with
  `MeshRegion[MeshOrderAlteration[mesh, 1]]`.

Timings referred to in the comments are for an
Intel(R) Core(TM) i7-7700 CPU @ 3.60 GHz (8 GB RAM)
with a `Benchmark[]` score of about 2.4.
