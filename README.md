# phd-code


Code for Conway's Ph.D., for Mathematica Version 12.

For the thesis, see [yawnoc/phd-thesis].

[yawnoc/phd-thesis]: https://github.com/yawnoc/phd-thesis


## BEFORE CLONING THIS REPOSITORY

**WARNING.**
**This repository contains upwards of 200 MB in objects**
(with a 500 MB working tree),
mostly compressed `` NDSolve`FEM` `` solutions,
which took a long time to compute.
In hindsight, it would have been better to store these files separately,
so as to not bloat this repository.

Including the large output files was a stupid way of ensuring
identical copies on a large number of computers.
While such a setup was reasonable using a USB stick as a shared Git remote,
fetching 200 MB per clone over the internet is *not* a good thing.

**To save bandwidth, CLONE [yawnoc/phd-code-no-output] INSTEAD,**
which is a debloated version of this repository with no output files.
Only download the large output files individually as needed.

[yawnoc/phd-code-no-output]: https://github.com/yawnoc/phd-code-no-output


## Figure fonts

LaTeX fonts in the figures were produced ***without*** Szabolcs's [MaTeX][].
While MaTeX is an excellent package outputting vector text labels,
one needs a TeX system and Ghostscript.
Instead I used a simple character mapping
along with Mathematica's own string formatting capabilities
(see the functions `LatinModernFont`, `LatinModernLabelStyle`,
`LatinModernFontStyle`, and `LaTeXStyle` in [`Conway.wl`](Conway.wl)).
This approach only requires the installation of the following fonts:

* "Latin Modern Math" `latinmodern-math.otf` <br>
  <http://www.gust.org.pl/projects/e-foundry/lm-math/download>

* "LM Roman 10" `lmroman10-*.otf` <br>
  <http://www.gust.org.pl/projects/e-foundry/latin-modern/download>

On Windows there are some [issues with spacing around Greek characters][greek]
due to a bug ... the solution is to use Linux instead.

Nevertheless, Linux appears to have a bug where
capital Greek letters fail to be italicised when exporting to PDF.
Therefore use Windows for [radiation-conduction-bvp.pdf],
[capillary-contact-angle.pdf], and [helmholtz-wedge-domain.pdf].

[greek]: https://mathematica.stackexchange.com/q/219745
[radiation-conduction-bvp.pdf]: radiation/general/radiation-conduction-bvp.pdf
[capillary-contact-angle.pdf]: capillary/general/capillary-contact-angle.pdf
[helmholtz-wedge-domain.pdf]: capillary/helmholtz/helmholtz-wedge-domain.pdf


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
Usually you want to run *only* the initialisation section.
This is the section I have called
"Initialisation section (always run this first)";
select its cell bracket and evaluate.
**Do NOT use "Evaluate Initialization Cells".**
(Cells in a `.wl` file are automatically styled as Code cells,
which are Initialization Cells by default.)

Timings referred to in the comments are for an
Intel(R) Core(TM) i7-7700 CPU @ 3.60 GHz (8 GB RAM)
with a `Benchmark[]` score of about 2.4.

[MaTeX]: https://github.com/szhorvat/MaTeX


## License

This repository is licensed under "MIT No Attribution" (MIT-0),
see [LICENSE](LICENSE).
