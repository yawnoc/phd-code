(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "conformal_triangle"};


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Conformal map (disk interior to triangle exterior)*)


(* ::Subsubsection:: *)
(*dz/d\[Zeta]*)


zOfZetaDerivative[zeta_, c_] := c / zeta^2 * (1 - zeta^3)^(2/3);
zOfZetaDerivative[\[FormalZeta], \[FormalCapitalC]]


(* ::Subsubsection:: *)
(*z(\[Zeta])*)


zOfZeta[zeta_, c_] := Integrate[zOfZetaDerivative[zeta, c], zeta] // Evaluate;
zOfZeta[\[FormalZeta], \[FormalCapitalC]]


(* ::Subsubsection:: *)
(*Multiplicative constant*)


constant = 1 / zOfZeta[1, 1];
{constant, constant // N}


zOfZetaDerivative[zeta_] := zOfZetaDerivative[zeta, constant] // Evaluate;
zOfZeta[zeta_] := zOfZeta[zeta, constant] // Evaluate;


(* Check *)
zOfZeta[1] == 1
zOfZeta'[\[FormalZeta]] == zOfZetaDerivative[\[FormalZeta]] // FullSimplify
