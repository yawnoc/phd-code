(* ::Package:: *)

(* ::Text:: *)
(*See r4 (manuscripts/radiation-4-polygon.pdf).*)


(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "polygon"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Schwarz--Christoffel transformation*)


(* ::Subsubsection:: *)
(*h_n(\[Zeta])*)


h[n_][zeta_] := Hypergeometric2F1[1/n, 2/n, 1 + 1/n, zeta ^ n];


(* ::Subsubsection:: *)
(*dz/d\[Zeta]*)


zMapDer[n_][zeta_] := 1 / (h[n][1] (1 - zeta ^ n) ^ (2/n)) // Evaluate;


(* ::Subsubsection:: *)
(*z = z(\[Zeta])*)


zMap[n_][zeta_] := zeta h[n][zeta] / h[n][1] // Evaluate;


(* ::Subsection:: *)
(*Italicised symbols*)


zIt = Italicised["z"];


(* ::Subsection:: *)
(*Global styles for plots*)


pointStyle = PointSize[Large];


(* ::Section:: *)
(*Schwarz--Christoffel transformation*)


(* ::Subsection:: *)
(*Forward transformation: disk (\[Zeta]-space) to polygon (z-space)*)


(* ::Text:: *)
(*See https://mathematica.stackexchange.com/a/159086*)


Table[
  Module[{rSamp = 4, phiSamp = 4},
    Show[
      ParametricPlot[r Exp[I phi] // zMap[n] // ReIm,
        {r, 0, 1}, {phi, 0, 2 Pi},
        Mesh -> {rSamp - 1, n phiSamp - 1},
        PlotPoints -> {rSamp, n phiSamp + 1},
        PlotOptions[Frame] // Evaluate
      ],
      Graphics @ {Directive[Red, pointStyle],
        Point @ CirclePoints[{1, 2 Pi / n}, n]
      }
    ]
  ] // Ex @ StringJoin["schwarz-christoffel-forward-", ToString[n], ".pdf"]
, {n, 3, 7}]
