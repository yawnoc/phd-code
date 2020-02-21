(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "cosine"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Known solution*)


(* ::Subsubsection:: *)
(*T*)


tKnown[b_][x_, y_] := 1 - b Cos[x] Cosh[y];


(* ::Subsubsection:: *)
(*P = \[PartialD]T/\[PartialD]x, Q = \[PartialD]T/\[PartialD]y*)


p[b_][x_, y_] := D[tKnown[b][x, y], x] // Evaluate;
q[b_][x_, y_] := D[tKnown[b][x, y], y] // Evaluate;


(* ::Subsection:: *)
(*Flux function F*)


f[a_, b_][x_, y_] := -tKnown[b][x, y]^4 / a // Evaluate;


(* ::Subsection:: *)
(*Viability \[CapitalPhi]*)


vi[a_, b_][x_, y_] := p[b][x, y]^2 + q[b][x, y]^2 - f[a, b][x, y]^2 // Evaluate;


(* ::Section:: *)
(*Algebra checks*)


(* ::Subsection:: *)
(*Check \[CapitalPhi]*)


With[{a = \[FormalCapitalA], b = \[FormalCapitalB], x = \[FormalX], y = \[FormalY]},
  vi[a, b][x, y] ==
    b^2 (Sin[x]^2 + Sinh[y]^2)
    - (1 - b Cos[x] Cosh[y])^8 / a^2
    // FullSimplify
]
