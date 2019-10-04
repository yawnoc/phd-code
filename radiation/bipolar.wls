#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*See r3 (manuscripts/radiation-3-bipolar.pdf).*)


(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< Curvilinear`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "bipolar"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Options for exported GIFS*)


gifOpts = Sequence[
  AnimationRepetitions -> Infinity,
  "DisplayDurations" -> 0.5
];


(* ::Subsection:: *)
(*Bipolar coordinates*)


(* ::Subsubsection:: *)
(*Effectively v = Infinity*)


vInfinity = 10;


(* ::Subsubsection:: *)
(*Contours*)


defContStyle = LightGray;


uContour[u_, style_: defContStyle] :=
  ParametricPlot[
    XYBipolar[u, v], {v, -vInfinity, vInfinity},
    Axes -> None,
    Exclusions -> If[Mod[u, 2 Pi] == 0, 0, Automatic],
    PlotStyle -> style
  ];


vContour[v_, style_: defContStyle] :=
  ParametricPlot[
    XYBipolar[u, v], {u, 0, 2 Pi},
    Axes -> None,
    PlotStyle -> style
  ];


uContours = uContour /@ Most[Range[0, 360, 30] Degree];
vContours = vContour /@ {-2, -1, -1/2, 0, 1/2, 1, 2};


(* ::Subsection:: *)
(*Global styles for plots*)


textStyle = Style[#, 18] &;


natStyle = Red;


nonStyle = Directive[Opacity[0.7], LightGray];
termStyle = Pink;
unphysStyle = Black;


upperStyle = Blue;
lowerStyle = Red;
convexStyle = Black;


uStyle = Blue;
vStyle = Red;


lighter = Directive[Opacity[0.5], #] &;
