#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*Wolfram Language scripts (.wls) are the way to go for Git.*)
(*https://mathematica.stackexchange.com/a/155268*)


SetDirectory @ NotebookDirectory[];


(* ::Section:: *)
(*Testing Conway.wl*)


<< Conway`


Plot[x^2, {x, 0, 1},
  PlotLabel -> BoxedLabel["This is a boxed label"],
  PlotOptions[Axes] // Evaluate
]


DefaultColours


Module[{y, iFun},
  iFun = NDSolveValue[{y'[x] == 0, y[0] == 1}, y, {x, 0, 1}, NoExtrapolation];
  {DomainStart[iFun], DomainEnd[iFun], iFun[1 + DomainEnd[iFun]]}
]


EmptyAxes[{0, 1}, {0, 1}]


EmptyFrame[{0, 1}, {0, 1}]


Module[{fun, x},
  fun = Cos[#] - # &;
  x = SeekRoot[fun, {0, 1}];
  {x, fun[x]}
]


SeekRootBisection[Cos[#] - # &, {0, 1}]
SeekRootBisection[Cos[#] - # &, {0, 1}, "ReturnIterations" -> True]


Way[0, 1]


Way["x", "y", "p"]
