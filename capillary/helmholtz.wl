(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< NDSolve`FEM`
<< Conway`
<< Curvilinear`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "helmholtz"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Wedge half-angle*)


(* ::Subsubsection:: *)
(*\[Alpha]*)


alpha = Pi/4;


(* ::Subsection:: *)
(*Known solution*)


(* ::Subsubsection:: *)
(*U*)


uKnown[x_, y_] := Exp[(-x+y)/Sqrt[2]] + Exp[(-x-y)/Sqrt[2]];


(* ::Subsubsection:: *)
(*P = \[PartialD]U/\[PartialD]x, Q = \[PartialD]U/\[PartialD]y*)


p = Derivative[1, 0] @ uKnown;
q = Derivative[0, 1] @ uKnown;


(* ::Subsection:: *)
(*Flux function F*)


f = 1;


(* ::Subsection:: *)
(*Viability \[CapitalPhi]*)


vi[x_, y_] := p[x, y]^2 + q[x, y]^2 - f^2 // Expand // Evaluate


(* ::Subsection:: *)
(*Terminal curve x = x(y)*)


xTerm[y_] := 1/Sqrt[2] * Log[2 Cosh[Sqrt[2] y]];


(* ::Subsection:: *)
(*Traced boundaries (arc-length parametrisation)*)


xDerivative[x_, y_, branchSign_: 1] :=
  Divide[
    -q[x, y] f + branchSign p[x, y] Re @ Sqrt @ vi[x, y],
    p[x, y]^2 + q[x, y]^2
  ];
yDerivative[x_, y_, branchSign_: 1] :=
  Divide[
    +p[x, y] f + branchSign q[x, y] Re @ Sqrt @ vi[x, y],
    p[x, y]^2 + q[x, y]^2
  ];


(* Arc-length parametrisation x == x(s), y == y(s). *)
xyTraced[
  {xInitial_, yInitial_}, sInitial_,
  {sStart_, sEnd_},
  sSign_: 1,
  branchSign_: 1,
  terminationVi_: 0
] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    NDSolveValue[
      {
        x'[s] == sSign * xDerivative[x[s], y[s], branchSign],
        y'[s] == sSign * yDerivative[x[s], y[s], branchSign],
        x[sInitial] == xInitial,
        y[sInitial] == yInitial,
        WhenEvent[vi[x[s], y[s]] < terminationVi, "StopIntegration"]
      }
      , {x, y}
      , {s, sStart, sEnd}
      , NoExtrapolation
    ]
  ];


(* ::Section:: *)
(*Algebra*)


(* ::Subsection:: *)
(*Check known solution*)


With[{x = \[FormalX], y = \[FormalY]},
  {
    (* Helmholtz equation *)
    Laplacian[uKnown[x, y], {x, y}] == uKnown[x, y],
    (* Boundary conditions *)
    XYPolar[1, alpha + Pi/2] . Grad[uKnown[x, y], {x, y}] == 1 /. {y -> x Tan[alpha]},
    XYPolar[1, -alpha - Pi/2] . Grad[uKnown[x, y], {x, y}] == 1 /. {y -> -x Tan[alpha]},
    (* Derivatives *)
    p[x, y] == -1/Sqrt[2] (Exp[(-x+y)/Sqrt[2]] + Exp[(-x-y)/Sqrt[2]]),
    q[x, y] ==  1/Sqrt[2] (Exp[(-x+y)/Sqrt[2]] - Exp[(-x-y)/Sqrt[2]]),
    Nothing
  } // FullSimplify
]


(* ::Subsection:: *)
(*Gradient squared*)


With[{x = \[FormalX], y = \[FormalY]},
  p[x, y]^2 + q[x, y]^2 // Expand
]


(* ::Subsection:: *)
(*Check viability*)


With[{x = \[FormalX], y = \[FormalY]},
  vi[x, y] == Exp[-Sqrt[2] (x-y)] + Exp[-Sqrt[2] (x+y)] - 1
    // FullSimplify
]


(* ::Subsection:: *)
(*Check terminal curve*)


With[{x = \[FormalX], y = \[FormalY]},
  vi[x, y] == 0 /. {x -> xTerm[y]}
    // FullSimplify
]


(* ::Section:: *)
(*Figure: wedge domain (helmholtz-wedge-domain)*)


Module[
  {
    wedgeXMax, wedgeYMax,
    clearanceProportion, xMin, xMax, yMax,
    alpha, angleMarkerRadius,
    normalVectorPosition, normalVectorLength,
    orthogonalityMarkerLength, orthogonalityMarkerStyle,
    wallThicknessCorrection,
    textStyle, testStyleOmega,
    dummyForTrailingCommas
  },
  (* Wedge dimensions *)
  wedgeXMax = 1;
  wedgeYMax = 1;
  (* Plot range *)
  clearanceProportion = 1/10;
  xMin = Way[0, wedgeXMax, -clearanceProportion];
  xMax = Way[0, wedgeXMax, 1 + clearanceProportion];
  yMax = Way[0, wedgeYMax, 1 + clearanceProportion];
  (* Wedge half-angle *)
  alpha = Pi/4;
  angleMarkerRadius = 1/5 wedgeXMax;
  (* Normal vector *)
  normalVectorPosition = 1/2 {wedgeXMax, wedgeYMax};
  normalVectorLength = Sqrt[2]/3 wedgeXMax;
  orthogonalityMarkerLength = 1/18 wedgeXMax;
  orthogonalityMarkerStyle = Directive[EdgeForm[Black], FaceForm[None]];
  wallThicknessCorrection = 0.1;
  (* Diagram *)
  textStyle = Style[#, 24] & @* LaTeXStyle;
  testStyleOmega = Style[#, 30] & @* LaTeXStyle;
  Show[
    EmptyAxes[{xMin, xMax}, {-yMax, yMax}
      , AspectRatio -> Automatic
      , AxesStyle -> Darker[Gray]
      , ImageSize -> 240
      , LabelStyle -> LatinModernLabelStyle[24]
      , Method -> {"AxesInFront" -> False}
      , Ticks -> None
    ],
    (* Wedge half-angle *)
    Graphics @ {
      Circle[{0, 0}, angleMarkerRadius, {0, alpha}]
    },
    Graphics @ {
      Text[
        "\[Alpha]" == SeparatedRow["/"]["\[Pi]", 4] // textStyle
        , XYPolar[angleMarkerRadius, alpha / 2]
        , {-1.2, -0.45}
      ]
    },
    (* Normal vector *)
    Graphics @ {Directive[GeneralStyle["Thick"], Arrowheads[Large]],
      Arrow @ {
        normalVectorPosition,
        normalVectorPosition + XYPolar[normalVectorLength, alpha + Pi/2]
      }
    },
    Graphics @ {
      Text[
        Embolden["n"] // textStyle
        , normalVectorPosition + XYPolar[normalVectorLength, alpha + Pi/2]
        , {0, -1}
      ]
    },
    Graphics @ {orthogonalityMarkerStyle,
      Rectangle[
        normalVectorPosition,
        normalVectorPosition + orthogonalityMarkerLength {1 - wallThicknessCorrection, 1}
      ] // Rotate[#, alpha, normalVectorPosition] &
    },
    (* Wedge domain *)
    Graphics @ {
      Text[
        "\[CapitalOmega]" // testStyleOmega
        , {3/4 wedgeXMax, 1/3 wedgeYMax}
        , {-1.2, -0.5}
      ]
    },
    (* Wedge domain boundary *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {{wedgeXMax, wedgeYMax}, {0, 0}, {wedgeXMax, -wedgeYMax}}
    },
    Graphics @ {
      Text[
        "\[PartialD]\[CapitalOmega]" // testStyleOmega
        , 4/5 {wedgeXMax, wedgeYMax}
        , {1, -1}
      ]
    },
    {}
  ]
] // Ex["helmholtz-wedge-domain.pdf"]


(* ::Section:: *)
(*Figure: known solution (helmholtz-wedge-solution)*)


Module[
  {
    wallHorizontalThickness, wallHeight,
    xMinWall, xMax, yMax, yMaxWall,
    verticalEdge, wedgeBaseVertices,
    solutionBoundingBoxVertices,
    makeBaseVertex3D,
    dummyForTrailingCommas
  },
  (* Wall properties *)
  wallHorizontalThickness = 0.25;
  wallHeight = 1.21 uKnown[0, 0];
  (* Plot range *)
  xMinWall = -wallHorizontalThickness;
  xMax = 4;
  yMax = xMax/2;
  yMaxWall = yMax + wallHorizontalThickness;
  (* Manually drawn wedge wall wedges *)
  verticalEdge @ {x_, y_} := Line @ {{x, y, 0}, {x, y, wallHeight}};
  wedgeBaseVertices = {
    {yMax - wallHorizontalThickness / 2, -yMax - wallHorizontalThickness / 2},
    {yMax, -yMax},
    {yMax, +yMax},
    {yMax - wallHorizontalThickness / 2, +yMax + wallHorizontalThickness / 2},
    Nothing
  };
  solutionBoundingBoxVertices = {
    {yMax, -yMax},
    {xMax, 0},
    {yMax, +yMax},
    Nothing
  };
  makeBaseVertex3D @ {x_, y_} := {x, y, 0};
  (* Plot *)
  Show[
    (* Known solution *)
    Plot3D[
      uKnown[x, y], {x, xMinWall, xMax}, {y, -yMaxWall, yMaxWall}
      , AxesEdge -> {{-1, -1}, {+1, -1}, {-1, -1}}
      , AxesLabel -> Italicise @ {"x", "y", "U"}
      , Boxed -> {Back, Bottom, Left}
      , BoxRatios -> Automatic
      , Filling -> 0
      , FillingStyle -> BoundaryTracingStyle["Solution3D"]
      , LabelStyle -> LatinModernLabelStyle[12]
      , Lighting -> GeneralStyle["AmbientLighting"]
      , PlotPoints -> 20
      , PlotStyle -> BoundaryTracingStyle["Solution3D"]
      , RegionFunction -> Function[{x, y},
          And[
            Abs[y] < x,
            Abs[y] < xMax - x
          ]
        ]
      , TicksStyle -> ConstantArray[LatinModernLabelStyle[8], 3]
      , ViewPoint -> {2.5, -1.1, 1.4}
    ],
    (* Wedge walls *)
    Plot3D[
      wallHeight, {x, xMinWall, xMax}, {y, -yMaxWall, yMaxWall}
      , Filling -> 0
      , FillingStyle -> BoundaryTracingStyle["Wall3D"]
      , Lighting -> GeneralStyle["AmbientLighting"]
      , Mesh -> None
      , PlotPoints -> 20
      , PlotStyle -> BoundaryTracingStyle["Wall3D"]
      , RegionFunction -> Function[{x, y},
          And[
            Abs[y] - wallHorizontalThickness < x < Abs[y],
            Abs[y] < xMax - x
          ]
        ]
    ],
    (* Manually drawn wedge wall edges *)
    Graphics3D @ {
      verticalEdge /@ wedgeBaseVertices,
      Line[makeBaseVertex3D /@ wedgeBaseVertices[[{1, 2}]]],
      Line[makeBaseVertex3D /@ wedgeBaseVertices[[{-1, -2}]]],
      Line @ {{0, 0, uKnown[0, 0]}, {0, 0, wallHeight}},
      {}
    },
    (* Manually drawn solution bounding box edges *)
    Graphics3D @ {
      Line[makeBaseVertex3D /@ solutionBoundingBoxVertices],
      Line @ {{xMax, 0, 0}, {xMax, 0, uKnown[xMax, 0]}},
      {}
    },
    {}
    , ImageSize -> 8 ImageSizeCentimetre
  ]
] // Ex["helmholtz-wedge-solution.png"
  , Background -> None
  , ImageResolution -> 4 BasicImageResolution
]


(* ::Section:: *)
(*Figure: non-viable domain (helmholtz-viable)*)


Module[
  {
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    dummyForTrailingCommas
  },
  (* Plot range *)
  xMax = 1;
  yMax = xMax;
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Plot *)
  Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}
      , ImageSize -> 240
      , LabelStyle -> LatinModernLabelStyle[12]
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {{xMaxMore, yMaxMore}, {0, 0}, {xMaxMore, -yMaxMore}}
    },
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0 && Abs[y] < x
      , {x, 0, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 5
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    {}
  ]
] // Ex["helmholtz-viable.pdf"]


(* ::Section:: *)
(*Figure: generic traced boundaries (helmholtz-traced-boundaries)*)


Module[
  {
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    xyStartList, xyTracedList,
    sStartPlotting,
    sEndPlotting,
    dummyForTrailingCommas
  },
  (* Plot range *)
  xMax = 1;
  yMax = xMax;
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Traced boundaries (upper branch) *)
  xyStartList = Table[XYPolar[r, -alpha], {r, Subdivide[rMax, 12]}];
  xyTracedList =
    Table[
      xyTraced[xyStart, 0, {0, 2 rMax}, -1]
      , {xyStart, xyStartList}
    ];
  (* Arc-length to start and end plotting at *)
  sStartPlotting[xy_] := DomainStart[xy];
  sEndPlotting[xy_] :=
    Module[{sStart, sEnd},
      sStart = DomainStart[xy];
      sEnd = DomainEnd[xy];
      If[sEnd < xMaxMore,
        sEnd,
        SeekRoot[
          xy[[1]][#] - xMaxMore &,
          {sStart, sEnd},
          5
        ]
      ]
    ];
  (* Plot *)
  Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}
      , ImageSize -> 240
      , LabelStyle -> LatinModernLabelStyle[12]
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {{xMaxMore, yMaxMore}, {0, 0}, {xMaxMore, -yMaxMore}}
    },
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0 && Abs[y] < x
      , {x, 0, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 5
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Traced boundaries *)
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // IncludeYReflection
          // Evaluate
        , {s, sStartPlotting[xy], sEndPlotting[xy]}
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {xy, xyTracedList}
    ],
    {}
  ]
] // Ex["helmholtz-traced-boundaries.pdf"]


(* ::Section:: *)
(*Figure: legend (helmholtz-legend)*)


Module[
  {
    legendLabelStyle,
    row,
    dummyForTrailingCommas
  },
  legendLabelStyle = LatinModernLabelStyle[14];
  row = Join[
    CurveLegend[
      BoundaryTracingStyle /@ {"Wall", "Terminal"},
      {"wedge wall", "terminal curve"}
      , LabelStyle -> legendLabelStyle
    ],
    RegionLegend[
      BoundaryTracingStyle /@ {"NonViable"},
      {"non\[Hyphen]viable domain"}
      , LabelStyle -> legendLabelStyle
    ],
    CurveLegend[
      BoundaryTracingStyle /@ {"Traced"},
      {"traced boundary"}
      , LabelStyle -> legendLabelStyle
    ],
    {}
  ];
  Grid[{row}
    , Alignment -> Left
    , Spacings -> {1, 0}
  ]
] // Ex["helmholtz-legend.pdf"]
