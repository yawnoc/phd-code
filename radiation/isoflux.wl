(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< FigureStyles`
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
tOfXY[x_, y_] := tOfZ[x + I y] // ComplexExpand // Evaluate;


(* ::Subsection:: *)
(*Slope dW/dz = dW/d\[Zeta] . d\[Zeta]/dz*)


wSlopeOfZ[z_] := wOfZeta'[zetaOfZ[z]] zetaOfZ'[z] // Evaluate;


(* ::Subsection:: *)
(*Flux function f*)


f = 1;


(* ::Subsection:: *)
(*Viability function \[CapitalPhi]*)


phiOfZ[z_] := Abs[wSlopeOfZ[z]]^2 - f // Evaluate;
phiOfZeta[zeta_] := phiOfZ[zOfZeta[zeta]] // Evaluate;


(* ::Subsection:: *)
(*Solution-positive region function*)


solutionPositiveRegionFunction =
  Function[{x, y}, tOfZ[x + I y] > 0 // Evaluate];


(* ::Subsection:: *)
(*Critical terminal point*)


zCritical = Module[{z}, z /. First @ Solve[phiOfZ[z] == 0 && z > 0, z, Reals]];
zetaCritical = zetaOfZ[zCritical];


(* ::Subsection:: *)
(*Traced boundary derivative d\[Zeta]/ds*)


zetaTracedDerivative[zeta_, sign_] :=
  Divide[
    I * f + (sign) Sqrt[phiOfZeta[zeta]],
    wOfZeta'[zeta]
  ];


(* ::Subsection:: *)
(*Traced boundary \[Zeta](s)*)


zetaTraced[zeta0_, {sStart_, sEnd_}, sign_, terminationPhi_: 0] :=
  Module[{zeta},
    NDSolveValue[
      {
        zeta'[s] == zetaTracedDerivative[zeta[s], sign],
        zeta[0] == zeta0,
        WhenEvent[phiOfZeta[zeta[s]] < 0, "StopIntegration"],
        WhenEvent[tOfZeta[zeta[s]] < 0, "StopIntegration"],
        Nothing
      }
      , zeta, {s, sStart, sEnd}
    ]
  ]


zetaTracedCriticalUpper = zetaTraced[zetaCritical, {-3, 3}, +1];
zetaTracedCriticalLower = zetaTraced[zetaCritical, {-3, 3}, -1];


(* ::Section:: *)
(*Checks*)


zOfZeta[\[FormalZeta]]
zetaOfZ[\[FormalZ]]


tOfZeta[\[FormalZeta]]
tOfZ[\[FormalZ]]


wSlopeOfZ[\[FormalZ]]


f


phiOfZ[\[FormalZ]]
phiOfZeta[\[FormalZeta]]


zCritical
zetaCritical


zetaTracedDerivative[\[FormalZeta], +1]
zetaTracedDerivative[\[FormalZeta], -1]


(* ::Section:: *)
(*Computer algebra*)


With[{s = \[FormalS], a0 = \[FormalA]},
  Module[
    {
      nMax, a, x, lhs, rhs,
      dummyForTrailingCommas
    },
    nMax = 8;
    a[0] = a0;
    x = Sum[a[n] s^n, {n, 0, nMax}] + O[s]^(nMax + 1);
    lhs = D[x, s];
    rhs = Sqrt[4 Sqrt[x^2 + s^2] - 1];
    Table[
      {
        n,
        a[n] =
          a[n] /. First @ Solve[
            SeriesCoefficient[lhs - rhs, {s, 0, n - 1}] == 0
              // FullSimplify[#, a0 > 0] &
              // Evaluate
            , a[n]
          ]
      }
      , {n, nMax}
    ]
  ]
] // TableForm


(* ::Section:: *)
(*Visualisation*)


(* ::Subsection:: *)
(*Known solution T*)


Plot3D[
  tOfZ[x + I y] // Evaluate
  , {x, 0, 3}
  , {y, -3, 3}
  , AspectRatio -> {1, 1, Automatic}
  , PlotRange -> {0, Automatic}
]


(* ::Subsection:: *)
(*Viability function \[CapitalPhi]*)


Plot3D[
  phiOfZ[x + I y] // Evaluate
  , {x, 0, 2}
  , {y, -2, 2}
  , AspectRatio -> {1, 1, Automatic}
  , PlotRange -> {0, Automatic}
  , RegionFunction -> solutionPositiveRegionFunction
]


(* ::Section:: *)
(*Boundary tracing*)


Show[
  EmptyFrame[
    {0, 1}, {-1, 1}
  ],
  (* Non-viable domain *)
  RegionPlot[
    phiOfZ[x + I y] < 0 && tOfZ[x + I y] > 0
    , {x, 0, 1}
    , {y, -1, 1}
    , BoundaryStyle -> None
    , PlotStyle -> BoundaryTracingStyle["NonViable"]
  ],
  (* Terminal curve *)
  ContourPlot[
    phiOfZ[x + I y] == 0
    , {x, 0, 1}
    , {y, -1, 1}
    , ContourStyle -> BoundaryTracingStyle["Terminal"]
    , RegionFunction -> solutionPositiveRegionFunction
  ],
  (* Solution contours *)
  ContourPlot[
    tOfZ[x + I y]
    , {x, 0, 1}
    , {y, -1, 1}
    , Contours ->
        Module[
          {
            tCriticalValue,
            divisionsUntoCritical,
            divisionsBeyondCritical,
            tSingleDivision,
            dummyForTrailingCommas
          },
          tCriticalValue = tOfZ[zCritical];
          divisionsUntoCritical = 2;
          divisionsBeyondCritical = 5;
          tSingleDivision = tCriticalValue / divisionsUntoCritical;
          tSingleDivision * Range[0, divisionsUntoCritical + divisionsBeyondCritical]
        ]
    , ContourShading -> None
    , ContourLabels -> None
  ],
  (* Traced boundaries (through critical terminal point) *)
  Table[
    ParametricPlot[
      zeta[s] // zOfZeta // ReIm // Evaluate
      , {s, DomainStart[zeta], DomainEnd[zeta]}
      , PlotStyle -> If[zeta === zetaTracedCriticalUpper, Blue, Red]
    ]
    , {zeta, {zetaTracedCriticalUpper, zetaTracedCriticalLower}}
  ],
  {}
]


(* ::Section:: *)
(*Sanity check*)


With[{x = \[FormalX], y = \[FormalY]},
  {
    {"Known solution T", tOfXY[x, y]},
    {"T along positive diagonal", tOfXY[x, x]},
    {"T along negative diagonal", tOfXY[x, -x]},
    {"Flux along vertical line", {1, 0}.Grad[tOfXY[x, y], {x, y}] /. {x -> Re[zCritical]}},
    Nothing
  }
] // TableForm
