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
