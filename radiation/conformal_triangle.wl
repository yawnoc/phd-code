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


(* ::Subsection:: *)
(*Known solution*)


(* ::Subsubsection:: *)
(*Analytical function, W(\[Zeta])*)


rho0 = 1/5;
w[zeta_] := Log[zeta / rho0] // Evaluate;


(* ::Subsubsection:: *)
(*Physical temperature, Re{W}*)


temperature[zeta_] := Re @ w[zeta] // Evaluate;


(* ::Subsection:: *)
(*Boundary condition*)


(* ::Subsubsection:: *)
(*Flux*)


a = 1;
flux[zeta_] := -temperature[zeta]^4 / a // Evaluate;


(* ::Subsubsection:: *)
(*Gradient squared*)


gradSquared[zeta_] := Abs[w'[zeta] / zOfZetaDerivative[zeta]] ^ 2 // Evaluate;
gradSquared[\[FormalZeta]] * constant^2


(* ::Subsubsection:: *)
(*Viability*)


viability[zeta_] := gradSquared[zeta] - flux[zeta]^2 // Evaluate;


(* ::Subsection:: *)
(*Hyperbolic critical terminal point*)


(* Principal point *)
angHyperbolic = Pi/3;
radHyperbolic = SeekRoot[
  viability[# Exp[I angHyperbolic]] &,
  {rho0, 1}
];
zetaHyperbolic = radHyperbolic * Exp[I angHyperbolic];


(* ::Subsection:: *)
(*Principal cubic root of unity*)


omega = Exp[I 2 Pi/3];
omega^3


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
  {xMin, xMax} = {yMin, yMax} = Norm[xyOfZeta[rho0]] {-1, 1};
  (* Contours (\[Zeta]-space) *)
  radValues = Subdivide[rho0, 1, 6];
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


(* ::Section:: *)
(*Visualise known solution*)


Module[
  {
    radMin, radMax,
    angMin, angMax,
    dummyForTrailingCommas
  },
  (* Plot range (\[Zeta]-space) *)
  {radMin, radMax} = {rho0, 1};
  {angMin, angMax} = {0, 2 Pi};
  (* Make plot *)
  ParametricPlot3D[
    Append[
      xyOfZeta[#],
      temperature[#]
    ] & [
      rad Exp[I ang]
    ]
    , {rad, radMin, radMax}
    , {ang, angMin, angMax}
    , BoxRatios -> {Automatic, Automatic, 2.5}
    , Exclusions -> None
  ]
]


(* ::Section:: *)
(*\[Zeta]-space visualisation*)


Module[
  {
    radValues, radMin, radMax,
    angValues, angMin, angMax,
    zetaContourStyle,
    dummyForTrailingCommas
  },
  (* Contours *)
  radValues = Subdivide[rho0, 1, 6];
  {radMin, radMax} = MinMax[radValues];
  angValues = Subdivide[0, 2 Pi, 12];
  {angMin, angMax} = MinMax[angValues];
  (* Styles *)
  zetaContourStyle = Lighter[Blue];
  (* Make plot *)
  Show[
    (* Radial contours *)
    ParametricPlot[
      Table[rad Exp[I ang] // ReIm, {rad, radValues}]
      , {ang, angMin, angMax}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> zetaContourStyle
    ],
    (* Azimuthal contours *)
    ParametricPlot[
      Table[rad Exp[I ang] // ReIm, {ang, angValues}]
      , {rad, radMin, radMax}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> zetaContourStyle
    ],
    (* Non-viable domain *)
    RegionPlot[
      viability[xx + I yy] < 0
      , {xx, -2 radMax, 2 radMax}
      , {yy, -2 radMax, 2 radMax}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 50
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Hyperbolic critical terminal points *)
    Graphics @ {
      Black, PointSize[Large],
      Point @ Table[zetaHyperbolic * omega^k // ReIm, {k, 0, 2}]
    },
    {}
  ]
]
