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


(* ::Subsubsection:: *)
(*z(\[Zeta]) with fixed constant*)


zOfZetaDerivative[zeta_] := zOfZetaDerivative[zeta, constant] // Evaluate;
zOfZeta[zeta_] := zOfZeta[zeta, constant] // Evaluate;


xyOfZeta[zeta_] := ReIm @ zOfZeta[zeta] // Evaluate;


(* Check *)
zOfZeta[1] == 1
zOfZeta'[\[FormalZeta]] == zOfZetaDerivative[\[FormalZeta]] // FullSimplify


(* ::Section:: *)
(*Visualise transformation*)


Module[
  {
    xMin, xMax, yMin, yMax,
    radValues, radMin, radMax,
    angValues, angMin, angMax,
    dummyForTrailingCommas
  },
  (* Plot range (z-space) *)
  {xMin, xMax} = {yMin, yMax} = 3 {-1, 1};
  (* Contours (\[Zeta]-space) *)
  radValues = Subdivide[0, 1, 8] // Rest;
  {radMin, radMax} = MinMax[radValues];
  angValues = Subdivide[0, 2 Pi, 12];
  {angMin, angMax} = MinMax[angValues];
  (* Make plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {yMin, yMax}
      , ImageSize -> 360
    ],
    (* Radial contours *)
    ParametricPlot[
      Table[rad Exp[I ang] // xyOfZeta, {rad, radValues}]
      , {ang, angMin, angMax}
      , PlotPoints -> 2
      , PlotRange -> Full
    ],
    (* Azimuthal contours *)
    ParametricPlot[
      Table[rad Exp[I ang] // xyOfZeta, {ang, angValues}]
      , {rad, radMin, radMax}
      , PlotPoints -> 2
      , PlotRange -> Full
    ],
    {}
  ]
]
