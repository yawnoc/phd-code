# phd-code

Code for Conway's Ph.D., for Mathematica Version&nbsp;12.

Licensed under "MIT No Attribution" (MIT-0), see [LICENSE](LICENSE).

See also: [phd-thesis](https://github.com/yawnoc/phd-thesis)

## Mathematica version

Everything should also work in Version&nbsp;11, except the following:

* Nonlinear PDEs / boundary conditions passed to `NDSolve` for
  `` NDSolve`FEM` ``.
  
* `Plot3D` where the domain is given in the form `Element[{x,y}, mesh]]`,
  where `mesh` is an `ElementMesh`.
  
  In Version&nbsp;11.0 this results in the error message
  > Plot3D::idomdim: {x,y} ∈ mesh does not have a valid dimension
    as a plotting domain.
  
  A [workaround] is to replace `mesh` with
  `MeshRegion[MeshOrderAlteration[mesh, 1]]`.

[workaround]: https://mathematica.stackexchange.com/a/180280

## Figure fonts

LaTeX fonts in the figures were produced ***without*** Szabolcs's [MaTeX][].
While MaTeX is an excellent package outputting vector text labels,
those labels cannot be selected and copied;
also one needs a TeX system and Ghostscript.
Instead I used a simple character mapping
along with Mathematica's own string formatting capabilities
(see the functions `LatinModernFont`, `LatinModernLabelStyle`,
`LatinModernFontStyle`, and `LaTeXStyle` in [`Conway.wl`](Conway.wl)).
This approach only requires the installation of the following fonts:

* "Latin Modern Math" `latinmodern-math.otf` <br>
  <http://www.gust.org.pl/projects/e-foundry/lm-math/download>

* "LM Roman 10" `lmroman10-regular.otf` <br>
  <http://www.gust.org.pl/projects/e-foundry/latin-modern/download>

## Workflow

Using Mathematica to do applied maths is an exploratory process.
Each section of the included notebooks contains explorations
which I thought useful at the time of writing.
Some are `Manipulate` visualisations,
some are static plots which I have exported as reminders for myself,
and some are complete dead ends.

**Don't simply run all the code in a notebook.**
This will not be very illuminating,
not to mention time-consuming since you will end up
doing computations which needn't be redone.
Usually you want to run *just* the initialisation section,
and after that, run a visualisation cell and then do some new exploration.

Timings referred to in the comments are for an
Intel(R) Core(TM) i7-7700 CPU @ 3.60 GHz (8 GB RAM)
with a `Benchmark[]` score of about 2.4.

[MaTeX]: https://github.com/szhorvat/MaTeX
