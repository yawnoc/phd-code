#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*Wolfram Language scripts (.wls) are the way to go for Git.*)
(*https://mathematica.stackexchange.com/a/155268*)


SetDirectory @ NotebookDirectory[]


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


"a + beta" // PrettyString["beta" -> "\[Beta]"]


SeekRootBisection[Cos[#] - # &, {0, 1}]
SeekRootBisection[Cos[#] - # &, {0, 1}, "ReturnIterations" -> True]


UniformRange[0, 1, 0.099]
UniformRange[0, 2 Pi, 5 Degree]


Way[0, 1]


Way["x", "y", "p"]


(* ::Section:: *)
(*Testing Curvilinear.wl*)


<< Curvilinear`


With[{r = \[FormalR], phi = \[FormalPhi], x = \[FormalX], y = \[FormalY]},
  {
    XPolar[r, phi],
    YPolar[r, phi],
    XYPolar[r, phi],
    RPolar[x, y],
    PhiPolar[x, y],
    RPhiPolar[x, y],
    HRPolar[r, phi],
    HPhiPolar[r, phi],
    ARPolar[r, phi],
    APhiPolar[r, phi],
    APhiPolar[r, phi] == Cross @ ARPolar[r, phi]
  } // TableForm[#, TableDepth -> 1] &
]


With[{u = \[FormalU], v = \[FormalV], x = \[FormalX], y = \[FormalY]},
  {
    XBipolar[u, v],
    YBipolar[u, v],
    XYBipolar[u, v],
    UBipolar[x, y],
    VBipolar[x, y],
    UVBipolar[x, y],
    HBipolar[u, v],
    HUBipolar[u, v],
    HVBipolar[u, v],
    AUBipolar[u, v],
    AVBipolar[u, v],
    AVBipolar[u, v] == Cross @ AUBipolar[u, v] // FullSimplify
  } // TableForm[#, TableDepth -> 1] &
]
