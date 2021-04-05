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
