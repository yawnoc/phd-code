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


(* ::Subsubsection:: *)
(*Principal cubic root of unity*)


omega = Exp[I 2 Pi/3];
omega^3


(* ::Subsubsection:: *)
(*Checks*)


zOfZeta[1] == 1
zOfZeta[omega^2] == omega
zOfZeta[omega] == omega^2
zOfZeta'[\[FormalZeta]] == zOfZetaDerivative[\[FormalZeta]] // FullSimplify


(* ::Subsection:: *)
(*Known solution*)


(* ::Subsubsection:: *)
(*Analytical function, W(\[Zeta])*)


rho0 = 0.2;
w[zeta_] := Log[zeta / rho0] // Evaluate;


(* ::Subsubsection:: *)
(*Physical temperature, Re{W}*)


temperature[zeta_] := Re @ w[zeta] // Evaluate;


(* ::Subsection:: *)
(*Boundary condition*)


(* ::Subsubsection:: *)
(*Flux*)


a = 1.5;
flux[zeta_] := -temperature[zeta] / a // Evaluate;


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
viabilityTolerance = 10^-6;
angHyperbolic = Pi/3;
radHyperbolic = SeekRootBisection[
  viability[# Exp[I angHyperbolic]] - viabilityTolerance &,
  {rho0, 1}
  , viabilityTolerance / 1000
];
zetaHyperbolic = radHyperbolic * Exp[I angHyperbolic];


viability[zetaHyperbolic]


(* ::Subsection:: *)
(*Traced boundaries*)


(* ::Subsubsection:: *)
(*Solvers for ODE*)


tracedRHS[zeta_, branchSign_: 1] :=
  Divide[
    I flux[zeta] + branchSign Sqrt[viability[zeta]],
    w'[zeta]
  ] // Evaluate;


sMax = 1.5 zOfZeta[rho0] // Ceiling;


tracedBoundary[zeta0_, branchSign_: 1, terminationViability_: 0] :=
  Module[{zeta},
    NDSolveValue[
      {
        zeta'[s] == tracedRHS[zeta[s], branchSign],
        zeta[0] == zeta0,
        WhenEvent[
          {
            Abs[zeta[s]] < rho0,
            Abs[zeta[s]] > 1,
            viability[zeta[s]] < terminationViability
          },
          "StopIntegration"
        ]
      },
      zeta, {s, -sMax, sMax}
      , NoExtrapolation
    ]
  ];


(* ::Subsubsection:: *)
(*Actual boundaries*)


(* ::Subsubsubsection:: *)
(*Starting points*)


tracedTypeList = {"general", "hyperbolic"};


zetaStartList["general"] =
  Module[{rad, angValues},
    rad = Way[rho0, radHyperbolic];
    angValues = 2 Pi / 3 {0.25, 0.5, 0.75};
    Table[rad Exp[I ang], {ang, angValues}]
  ];


zetaStartList["hyperbolic"] = {zetaHyperbolic};


(* ::Subsubsubsection:: *)
(*Boundaries*)


Table[
  zetaTracedList[type, branchSign] =
    Table[tracedBoundary[zeta0, branchSign], {zeta0, zetaStartList[type]}
    ];
  , {type, tracedTypeList}
  , {branchSign, {-1, 1}}
];


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
    upperStyle, lowerStyle, hyperbolicStyle,
    dummyForTrailingCommas
  },
  (* Contours *)
  radValues = Subdivide[rho0, 1, 6];
  {radMin, radMax} = MinMax[radValues];
  angValues = Subdivide[0, 2 Pi, 12];
  {angMin, angMax} = MinMax[angValues];
  (* Styles *)
  zetaContourStyle = Nest[Lighter, Blue, 3];
  hyperbolicStyle = Directive[Thick, Black];
  upperStyle = Blue;
  lowerStyle = Red;
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
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          zeta[s] * omega^k // ReIm // Evaluate
          , {s, DomainStart[zeta], DomainEnd[zeta]}
          , PlotStyle -> Which[
              type == "hyperbolic", hyperbolicStyle,
              branchSign == 1, upperStyle,
              branchSign == -1, lowerStyle,
              True, Automatic
            ]
        ]
        , {zeta, zetaTracedList[type, branchSign]}
      ]
      , {type, tracedTypeList}
      , {branchSign, {-1, 1}}
      , {k, 0, 2}
    ],
    {}
  ]
]


(* ::Section:: *)
(*z-space visualisation*)


Module[
  {
    xMin, xMax, yMin, yMax,
    radValues, radMin, radMax,
    angValues, angMin, angMax,
    zetaContourStyle,
    upperStyle, lowerStyle, hyperbolicStyle,
    dummyForTrailingCommas
  },
  (* Plot range (z-space) *)
  {xMin, xMax} = {yMin, yMax} = Norm[xyOfZeta[rho0]] {-1, 1};
  (* Contours (\[Zeta]-space) *)
  radValues = Subdivide[rho0, 1, 6];
  {radMin, radMax} = MinMax[radValues];
  angValues = Subdivide[0, 2 Pi, 12];
  {angMin, angMax} = MinMax[angValues];
  (* Styles *)
  zetaContourStyle = Nest[Lighter, Blue, 3];
  hyperbolicStyle = Directive[Thick, Black];
  upperStyle = Blue;
  lowerStyle = Red;
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
      , PlotStyle -> zetaContourStyle
    ],
    (* Azimuthal contours *)
    ParametricPlot[
      Table[rad Exp[I ang] // xyOfZeta, {ang, angValues}]
      , {rad, radMin, radMax}
      , PlotPoints -> 2
      , PlotRange -> Full
      , PlotStyle -> zetaContourStyle
    ],
    (* Hyperbolic critical terminal points *)
    Graphics @ {
      Black, PointSize[Large],
      Point @ Table[zetaHyperbolic * omega^k // xyOfZeta, {k, 0, 2}]
    },
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          zeta[s] * omega^k // xyOfZeta // Evaluate
          , {s, DomainStart[zeta], DomainEnd[zeta]}
          , PlotPoints -> 2
          , PlotStyle -> Which[
              type == "hyperbolic", hyperbolicStyle,
              branchSign == 1, upperStyle,
              branchSign == -1, lowerStyle,
              True, Automatic
            ]
        ]
        , {zeta, zetaTracedList[type, branchSign]}
      ]
      , {type, tracedTypeList}
      , {branchSign, {-1, 1}}
      , {k, 0, 2}
    ],
    {}
  ]
]
