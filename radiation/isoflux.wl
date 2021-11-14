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


zOfZeta[zeta_] := zeta^(1/2);
zetaOfZ[z_] := z^2;


(* ::Subsection:: *)
(*Analytic function W*)


wOfZeta[zeta_] := zeta;


(* ::Subsection:: *)
(*Known solution T = Re[W]*)


tOfZeta[zeta_] := Re[wOfZeta[zeta]] // Evaluate;
tOfZ[z_] := tOfZeta[zetaOfZ[z]] // Evaluate;


(* ::Section:: *)
(*Visualisation*)


(* ::Subsection:: *)
(*Known solution*)


Plot3D[
  tOfZ[x + I y] // Evaluate
  , {x, 0, 3}
  , {y, -3, 3}
  , AspectRatio -> {1, 1, Automatic}
  , PlotRange -> {0, Automatic}
]
