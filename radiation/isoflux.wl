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
(*Terminal point of given phase*)


zTerminal[phase_] := Module[{z}, z /. First @ Solve[phiOfZ[z] == 0 && Arg[z] == phase, z, Complexes]];
zetaTerminal[phase_] := zetaOfZ @ zTerminal[phase];


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
        WhenEvent[phiOfZeta[zeta[s]] < terminationPhi, "StopIntegration"],
        WhenEvent[tOfZeta[zeta[s]] < 0, "StopIntegration"],
        Nothing
      }
      , zeta, {s, sStart, sEnd}
    ]
  ]


(* ::Subsubsection:: *)
(*(from walls)*)


yWallList = Join[
  Subdivide[1/2, 1, 5],
  Subdivide[1/2, 1, 5] // Minus // Reverse,
  {}
];


zetaTracedWallListUpper =
  Module[{zeta0},
    Table[
      zeta0 = zetaOfZ[Abs[y] + I y];
      zetaTraced[zeta0, {0, 3}, +1]
      , {y, yWallList // Most}
    ]
  ];


zetaTracedWallListLower =
  Module[{zeta0},
    Table[
      zeta0 = zetaOfZ[Abs[y] + I y];
      zetaTraced[zeta0, {0, -3}, -1]
      , {y, yWallList // Rest}
    ]
  ];


(* ::Subsubsection:: *)
(*(from terminal curve)*)


phiTracedTerminalList = Pi/4 {-1, -3/4, -1/2, 1/2, 3/4, 1};


zetaTracedTerminalListUpper =
  Module[{zeta0},
    Table[
      zeta0 = zetaTerminal[phase];
      zetaTraced[zeta0, Sign[phase] {0, 3}, +1]
      , {phase, phiTracedTerminalList // Rest}
    ]
  ];


zetaTracedTerminalListLower =
  Module[{zeta0},
    Table[
      zeta0 = zetaTerminal[phase];
      zetaTraced[zeta0, Sign[phase] {0, 3}, -1]
      , {phase, phiTracedTerminalList // Most}
    ]
  ];


(* ::Subsubsection:: *)
(*(from axis)*)


zetaAxisPoint = zetaOfZ[7/10];


zetaTracedAxisUpper = zetaTraced[zetaAxisPoint, {-3, 3}, +1];
zetaTracedAxisLower = zetaTraced[zetaAxisPoint, {-3, 3}, -1];


(* ::Subsubsection:: *)
(*(through critical terminal points)*)


zetaTracedCriticalUpper = zetaTraced[zetaCritical, {-3, 3}, +1];
zetaTracedCriticalLower = zetaTraced[zetaCritical, {-3, 3}, -1];


(* ::Subsection:: *)
(*Styles*)


upperStyle = Blue;
lowerStyle = Red;


criticalStyle = Directive[Yellow, PointSize[Large]];


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


zTerminal[0] == zCritical


zetaTracedDerivative[\[FormalZeta], +1]
zetaTracedDerivative[\[FormalZeta], -1]


(* ::Section:: *)
(*Computer algebra*)


(* ::Subsection:: *)
(*Initial \[Zeta] != 1/4*)


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
        Subscript[Italicise["a"], n],
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
  ] // ReplaceAll[a0 -> Subscript[Italicise["a"], 0]]
] // TableForm // Ex["isoflux-ordinary-series.pdf"]


(* ::Subsection:: *)
(*Initial \[Zeta] = 1/4 and s > 0*)


With[{s = \[FormalS], a0 = \[FormalA]},
  Module[
    {
      nMax, a, x, lhs, rhs,
      dummyForTrailingCommas
    },
    nMax = 16;
    a[0] = a0;
    x = 1/4 + 2 s^2 + Sum[a[n] s^n, {n, 4, nMax, 2}] + O[s]^(nMax + 2);
    lhs = D[x, s];
    rhs = Sqrt[4 Sqrt[x^2 + s^2] - 1];
    Table[
      {
        Subscript[Italicise["a"], n],
        a[n] =
          a[n] /. First @ Solve[
            SeriesCoefficient[lhs - rhs, {s, 0, n - 1}] == 0
              // FullSimplify[#, s > 0] &
              // Evaluate
            , a[n]
          ]
      }
      , {n, 4, nMax, 2}
    ]
  ] // ReplaceAll[a0 -> Subscript[Italicise["a"], 0]]
] // TableForm // Ex["isoflux-special-series-positive.pdf"]


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


(* ::Section:: *)
(*Figure: isoflux BVP (isoflux-bvp)*)


Module[
  {
    yZeroValueUpper, yZeroValueLower,
    xRight, yTopRight, yBottomRight,
    xUnitFlux,
    textStyle,
    dummyForTrailingCommas
  },
  (* Zero-value boundaries *)
  yZeroValueUpper[x_] := x;
  yZeroValueLower[x_] := -yZeroValueUpper[x];
  (* Corner coordinates *)
  xRight = 1;
  yTopRight = yZeroValueUpper[xRight];
  yBottomRight = yZeroValueLower[xRight];
  (* Unit-flux boundary *)
  xUnitFlux[y_] := xRight + 1/7 Sin[Pi * y];
  (* Plot *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  Show[
    Graphics[],
    (* Zero-value boundaries *)
    Plot[
      {yZeroValueUpper[x], yZeroValueLower[x]}
      , {x, 0, xRight}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Contour"]
    ],
    Graphics @ {
      Text[
        "zero\[Hyphen]value" // textStyle
        , {#, yZeroValueUpper[#]}
        , {0.05, -1.1}
        , {1, yZeroValueUpper'[#]}
      ] & [xRight/2],
      Text[
        "zero\[Hyphen]value" // textStyle
        , {#, yZeroValueLower[#]}
        , {0.05, 0.75}
        , {1, yZeroValueLower'[#]}
      ] & [xRight/2],
      {}
    },
    (* Unit-flux boundary *)
    ParametricPlot[
      {xUnitFlux[y], y}
      , {y, yBottomRight, yTopRight}
      , PlotPoints -> 3
      , PlotStyle -> BoundaryTracingStyle["Traced"]
    ],
    Graphics @ {
      Text[
        "unit\[Hyphen]flux" // textStyle
        , {xUnitFlux[#], #}
        , {0, 0.9}
        , {xUnitFlux'[#], 1}
      ] & [0.03],
      {}
    },
    {}
    , ImageSize -> 0.25 ImageSizeTextWidth
  ]
] // Ex["isoflux-bvp.pdf"]


(* ::Section:: *)
(*Figure: grid (isoflux-grid-*-space)*)


(* ::Subsection:: *)
(*\[Zeta]-space*)


Module[
  {
    xxValues,
    xxMin, xxMax, yyMax,
    axesLabel,
    dummyForTrailingCommas
  },
  (* Grid lines (Re[zeta] == const) *)
  xxValues = Subdivide[0, 3/2, 8] // Rest;
  (* Make plot *)
  xxMin = 0;
  xxMax = 1.2;
  yyMax = 1;
  axesLabel[string_] := Row @ {string, " ", "\[Zeta]" // LaTeXStyle};
  Show[
    EmptyFrame[
      {xxMin, xxMax}, {-yyMax, yyMax}
      , FrameLabel -> {
          axesLabel["Re"] // Margined @ {{0, 0}, {0, -10}},
          axesLabel["Im"] // Margined @ {{0, -4}, {0, 0}}
        }
      , FrameLabel -> {
          axesLabel["Re"] // Margined @ {{0, 0}, {0, -10}},
          axesLabel["Im"] // Margined @ {{0, -4}, {0, 0}}
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle[LabelSize["Axis"] - 1]
      , PlotRangePadding -> {Scaled[0.05], Scaled[0.03]}
    ],
    (* Generic grid lines *)
    Graphics @ {
      BoundaryTracingStyle["Background"],
      Table[
        Line @ {{xx, -#}, {xx, #}} & [1.1 yyMax]
      , {xx, xxValues}]
    },
    (* Zero-value boundary *)
    Graphics @ {
      BoundaryTracingStyle["Contour"],
      Table[
        Line @ {{xx, -#}, {xx, #}} & [1.1 yyMax]
      , {xx, {0}}]
    },
    {}
    , ImageSize -> 0.4 ImageSizeTextWidth
  ]
] // Ex["isoflux-grid-zeta-space.pdf"]


(* ::Subsection:: *)
(*z-space*)


Module[
  {
    xxValues, xValues, xCurve,
    xMin, xMax, yMax,
    axesLabel,
    dummyForTrailingCommas
  },
  (* Grid lines (Re[zeta] == const) *)
  (*
      x_0^2 == x^2 - y^2
      x == sqrt(x_0^2 + y^2)
  *)
  xxValues = Subdivide[0, 3/2, 8] // Rest;
  xValues = Sqrt[xxValues];
  xCurve[x0_][y_] := Sqrt[x0^2 + y^2];
  (* Make plot *)
  xMin = 0;
  xMax = 1.2;
  yMax = 1;
  axesLabel[string_] := Row @ {string, "\[ThinSpace]", Italicise["z"]};
  Show[
    EmptyFrame[
      {xMin, xMax}, {-yMax, yMax}
      , FrameLabel -> {
          axesLabel["Re"] // Margined @ {{0, 0}, {0, -15}},
          axesLabel["Im"] // Margined @ {{0, -2}, {0, 0}}
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle[LabelSize["Axis"] - 1]
    ],
    (* Generic grid lines (curves) *)
    ParametricPlot[
      Table[{xCurve[x0][y], y}, {x0, xValues}] // Evaluate
      , {y, -1.1 yMax, 1.1 yMax}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Background"]
    ],
    (* Zero-value boundary *)
    ParametricPlot[
      Table[{xCurve[x0][y], y}, {x0, {0}}] // Evaluate
      , {y, -1.1 yMax, 1.1 yMax}
      , Exclusions -> None
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Contour"]
    ],
    {}
    , ImageSize -> 0.4 ImageSizeTextWidth
  ]
] // Ex["isoflux-grid-z-space.pdf"]


(* ::Section:: *)
(*Figure: traced boundaries (isoflux-traced-boundaries-*-space)*)


(* ::Subsection:: *)
(*\[Zeta]-space*)


Module[
  {
    xxMin, xxMax, yyMax,
    axesLabel,
    dummyForTrailingCommas
  },
  xxMin = 0;
  xxMax = 1.2;
  yyMax = 1;
  axesLabel[string_] := Row @ {string, " ", "\[Zeta]" // LaTeXStyle};
  Show[
    EmptyFrame[
      {xxMin, xxMax}, {-yyMax, yyMax}
      , FrameLabel -> {
          axesLabel["Re"] // Margined @ {{0, 0}, {0, -10}},
          axesLabel["Im"] // Margined @ {{0, -4}, {0, 0}}
        }
      , FrameLabel -> {
          axesLabel["Re"] // Margined @ {{0, 0}, {0, -10}},
          axesLabel["Im"] // Margined @ {{0, -4}, {0, 0}}
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle[LabelSize["Axis"] - 1]
      , PlotRangePadding -> Scaled[0.03]
    ],
    (* Non-viable domain *)
    RegionPlot[
      phiOfZeta[xx + I yy] < 0 && tOfZeta[xx + I yy] > 0
      , {xx, xxMin, xxMax}
      , {yy, -yyMax, yyMax}
      , BoundaryStyle -> None
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Traced boundaries (from walls) *)
    Table[
      ParametricPlot[
        zeta[s] // ReIm // Evaluate
        , {s, DomainStart[zeta], DomainEnd[zeta]}
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {zeta, zetaTracedWallListUpper}
    ],
    Table[
      ParametricPlot[
        zeta[s] // ReIm // Evaluate
        , {s, DomainStart[zeta], DomainEnd[zeta]}
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {zeta, zetaTracedWallListLower}
    ],
    (* Traced boundaries (from terminal curve) *)
    Table[
      ParametricPlot[
        zeta[s] // ReIm // Evaluate
        , {s, DomainStart[zeta], DomainEnd[zeta]}
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {zeta, zetaTracedTerminalListUpper}
    ],
    Table[
      ParametricPlot[
        zeta[s] // ReIm // Evaluate
        , {s, DomainStart[zeta], DomainEnd[zeta]}
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {zeta, zetaTracedTerminalListLower}
    ],
    (* Traced boundaries (from axis) *)
    Table[
      ParametricPlot[
        zeta[s] // ReIm // Evaluate
        , {s, DomainStart[zeta], DomainEnd[zeta]}
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {zeta, {zetaTracedAxisUpper, zetaTracedAxisLower}}
    ],
    (* Traced boundaries (through critical terminal point) *)
    Table[
      ParametricPlot[
        zeta[s] // ReIm // Evaluate
        , {s, DomainStart[zeta], DomainEnd[zeta]}
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {zeta, {zetaTracedCriticalUpper, zetaTracedCriticalLower}}
    ],
    (* Zero-value boundary *)
    Graphics @ {
      BoundaryTracingStyle["Contour"],
      Line @ {{0, -#}, {0, #}} & [1.1 yyMax]
    },
    {}
    , ImageSize -> 0.4 ImageSizeTextWidth
  ]
] // Ex["isoflux-traced-boundaries-zeta-space.pdf"]


(* ::Subsection:: *)
(*z-space*)


Module[
  {
    xMin, xMax, yMax,
    axesLabel,
    dummyForTrailingCommas
  },
  xMin = 0;
  xMax = 1.2;
  yMax = 1;
  axesLabel[string_] := Row @ {string, "\[ThinSpace]", Italicise["z"]};
  Show[
    EmptyFrame[
      {xMin, xMax}, {-yMax, yMax}
      , FrameLabel -> {
          axesLabel["Re"] // Margined @ {{0, 0}, {0, -15}},
          axesLabel["Im"] // Margined @ {{0, -2}, {0, 0}}
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle[LabelSize["Axis"] - 1]
    ],
    (* Non-viable domain *)
    RegionPlot[
      phiOfZ[x + I y] < 0 && tOfZ[x + I y] > 0
      , {x, xMin, xMax}
      , {y, -yMax, yMax}
      , BoundaryStyle -> None
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Traced boundaries (from walls) *)
    Table[
      ParametricPlot[
        zeta[s] // zOfZeta // ReIm // Evaluate
        , {s, DomainStart[zeta], DomainEnd[zeta]}
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {zeta, zetaTracedWallListUpper}
    ],
    Table[
      ParametricPlot[
        zeta[s] // zOfZeta // ReIm // Evaluate
        , {s, DomainStart[zeta], DomainEnd[zeta]}
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {zeta, zetaTracedWallListLower}
    ],
    (* Traced boundaries (from terminal curve) *)
    Table[
      ParametricPlot[
        zeta[s] // zOfZeta // ReIm // Evaluate
        , {s, DomainStart[zeta], DomainEnd[zeta]}
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {zeta, zetaTracedTerminalListUpper}
    ],
    Table[
      ParametricPlot[
        zeta[s] // zOfZeta // ReIm // Evaluate
        , {s, DomainStart[zeta], DomainEnd[zeta]}
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {zeta, zetaTracedTerminalListLower}
    ],
    (* Traced boundaries (from axis) *)
    Table[
      ParametricPlot[
        zeta[s] // zOfZeta // ReIm // Evaluate
        , {s, DomainStart[zeta], DomainEnd[zeta]}
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {zeta, {zetaTracedAxisUpper, zetaTracedAxisLower}}
    ],
    (* Traced boundaries (through critical terminal point) *)
    Table[
      ParametricPlot[
        zeta[s] // zOfZeta // ReIm // Evaluate
        , {s, DomainStart[zeta], DomainEnd[zeta]}
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {zeta, {zetaTracedCriticalUpper, zetaTracedCriticalLower}}
    ],
    (* Zero-value boundary *)
    Graphics @ {
      BoundaryTracingStyle["Contour"],
      Line @ {{#, #}, {0, 0}, {#, -#}} & [1.1 xMax]
    },
    {}
    , ImageSize -> 0.4 ImageSizeTextWidth
  ]
] // Ex["isoflux-traced-boundaries-z-space.pdf"]


(* ::Subsection:: *)
(*Legend*)


Module[{legendLabelStyle},
  legendLabelStyle = LatinModernLabelStyle @ LabelSize["Legend"];
  GraphicsGrid[
    Transpose @ {
      CurveLegend[
        {BoundaryTracingStyle["Traced"]},
        {"traced boundary"}
        , LabelStyle -> legendLabelStyle
      ],
      RegionLegend[
        {BoundaryTracingStyle["NonViable"]},
        {"non\[Hyphen]viable domain"}
        , LabelStyle -> legendLabelStyle
      ],
      Nothing
    }
    , Alignment -> Left
    , ImageSize -> 0.65 ImageSizeTextWidth
    , ItemAspectRatio -> 0.15
  ]
] // Ex["isoflux-traced-boundaries-legend.pdf"]
