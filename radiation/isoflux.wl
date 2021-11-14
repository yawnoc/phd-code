(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "isoflux"};


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Conformal mapping*)


zMap[zeta_] := zeta^(1/2);
zetaMap[z_] := z^2;


(* ::Subsection:: *)
(*Analytic function W*)


w[zeta_] := zeta;


(* ::Subsection:: *)
(*Known solution T = Re[W]*)


t[z_] := Re[w[zetaMap[z]]] // Evaluate;


t[\[FormalZ]]


(* ::Section:: *)
(*Visualisation*)


(* ::Subsection:: *)
(*Known solution*)


Plot3D[
  t[x + I y] // Evaluate
  , {x, 0, 3}
  , {y, -3, 3}
  , AspectRatio -> {1, 1, Automatic}
  , PlotRange -> {0, Automatic}
]
