(* ::Package:: *)

(* ::Text:: *)
(*See r1 (manuscripts/radiation-1-plane.pdf).*)


(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "plane"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Traced boundary y = y(x)*)


(* ::Text:: *)
(*See (r1.17) to (r1.20) (Page r1-2 of manuscripts/radiation-1-plane.pdf).*)


(* NOTE: making yTraDer negative was a really dumb decision. Too late to change it though. *)
yTraDer = Function[{x}, -x^4 / Sqrt[1 - x^8]];
yTra = Function[{x}, Integrate[yTraDer[x], x] // Evaluate];
yTraMax = Abs @ yTra[1];


With[{x = \[FormalX]},
  {yTraDer[x], yTra[x], yTraMax}
]


(* ::Subsection:: *)
(*Self-incident boundary ratio*)


(* ::Subsubsection:: *)
(*Exact expression*)


(* ::Text:: *)
(*See (r6.15) (Page r6-2 of manuscripts/radiation-6-self.pdf).*)


boundaryRatioIntegrand[x_, xx_] :=
  Module[{y, yDer},
    y = -yTra[#] &;
    yDer = -yTraDer[#] &;
    Divide[
      Times[
        xx^4,
        -(x - xx) yDer[xx] + (y[x] - y[xx]),
        +(x - xx) yDer[x]  - (y[x] - y[xx])
      ],
      Times[
        2,
        ((x - xx)^2 + (y[x] - y[xx])^2) ^ (3/2),
        Sqrt[1 + yDer[x]^2]
      ]
    ]
  ] // Evaluate;


boundaryRatio[x1_, x2_][x_] :=
  1/x^4 * NIntegrate[boundaryRatioIntegrand[x, xx], {xx, x1, x2}];


(* ::Subsubsection:: *)
(*Crude upper bound*)


(* ::Text:: *)
(*See (r6.22) (Page r6-3 of manuscripts/radiation-6-self.pdf).*)


boundaryRatioBoundIntegrand[x_, xx_] :=
  Module[{y},
    y = -yTra[#] &;
    Divide[
      xx^4 * (x - xx)^4,
      ((x - xx)^2 + (y[x] - y[xx])^2) ^ (3/2)
    ]
  ] // Evaluate;


boundaryRatioBoundPrefactor[x1_, x2_][x_] :=
  Module[{y, yDer, yDer2},
    y = -yTra[#] &;
    yDer = -yTraDer[#] &;
    yDer2 = yDer';
    Divide[
      (yDer2 @ Max[x1, x2, x])^2,
      8 x^4 * Sqrt[1 + yDer[x]^2]
    ]
  ] // Evaluate;


boundaryRatioBound[x1_, x2_][x_] :=
  Times[
    boundaryRatioBoundPrefactor[x1, x2][x],
    NIntegrate[boundaryRatioBoundIntegrand[x, xx], {xx, x1, x2}]
  ]


(* ::Subsubsection:: *)
(*Ultra-crude upper bound*)


(* ::Text:: *)
(*See (r6.26) (Page r6-4 of manuscripts/radiation-6-self.pdf).*)


boundaryRatioBoundUltra[x1_, x2_] /; x1 < x2 :=
  Module[{y, yDer, yDer2},
    y = -yTra[#] &;
    yDer = -yTraDer[#] &;
    yDer2 = yDer';
    Divide[
      (x2 - x1)^2 * (x2/x1)^4 * yDer2[x2]^2,
      8 (1 + yDer[x1]^2) ^ 2
    ] // N
  ] // Evaluate;


(* ::Subsection:: *)
(*Parameters for figures*)


(* Avoid confusion (note yTra == -2F1(x; ...)) *)
yTraUpper[c_][x_] := c + yTra[x];
yTraLower[c_][x_] := c - yTra[x];


(* ::Subsubsection:: *)
(*Patched portions of traced boundaries*)


(*
  REMEMBER:
    the positive signed branch is the lower;
    the negative signed branch is the upper.
  This is because the sign in the tracing equation is \[MinusPlus].
  Therefore in a spike (shaped like <),
  the lower-branch curve is actually higher.
  However, writing each curve as y == c \[MinusPlus] 2F1(x; ...),
  note that the lower-branch curve has lower c.
  Now, to generate assorted patchings of traced boundaries,
  we choose the corners (x_1, y_1), ..., (x_n, y_n)
  with y_1 < ... < y_n
  and determine the intersections
  between lower-branch(i) and upper-branch(i+1)
  for i == 1, ..., n - 1.
  The corners (x_i, y_i) cannot be chosen arbitrarily,
  so don't be dumb.
 *)


(* Lists of corners *)
patchedCornerList = Association[
  "regular-long" -> Table[{0.4, y}, {y, Subdivide[-0.6, 0.6, 3]}],
  "regular-short" -> Table[{0.8, y}, {y, Subdivide[-0.7, 0.7, 8]}],
  "irregular" -> {
    {0.6, -0.54},
    {0.81, -0.435},
    {0.77, -0.15},
    {0.23, -0.06},
    {0.45, 0.},
    {0.7, 0.15},
    {0.3, 0.6}
  }
];


patchedIdList = Keys[patchedCornerList];


Table[
  {patchedCornerXList[id], patchedCornerYList} =
    Transpose @ patchedCornerList[id]
, {id, patchedIdList}];


(*
  Lists of constants {cUpper, cLower}
  for the pair of traced boundaries through each corner point
 *)
Table[
  Module[{xCorner, yCorner, cUpper, cLower},
    {patchedCUpperList[id], patchedCLowerList[id]} =
      Transpose @ Table[
        {xCorner, yCorner} = corner;
        cUpper = yCorner - yTra[xCorner];
        cLower = yCorner + yTra[xCorner];
        {cUpper, cLower}
      , {corner, patchedCornerList[id]}]
  ]
, {id, patchedIdList}];


(*
  Lists of intersections between lower-branch(i) and upper-branch(i+1)
 *)
Table[
  patchedCornerNum[id] = Length @ patchedCornerList[id];
  Module[{cUpper, cLower, xInt, yInt},
    {patchedIntXList[id], patchedIntYList[id]} =
      Transpose @ Table[
        cLower = patchedCLowerList[id][[i]];
        cUpper = patchedCUpperList[id][[i + 1]];
        xInt = SeekRoot[
          yTraUpper[cUpper][#] - yTraLower[cLower][#] &,
          {0, 1}
        ];
        yInt = yTraUpper[cUpper][xInt];
        {xInt, yInt}
      , {i, patchedCornerNum[id] - 1}]
  ]
, {id, patchedIdList}];


(* ::Subsubsection:: *)
(*Corners for constructing domains*)


(*
  Corners used to build domains (specified by index i)
  among each list of corners
*)
domainCornerRangeList["regular-long"] = {{2, 3}};
domainCornerRangeList["regular-short"] = {{3, 7}};
domainCornerRangeList["irregular"] = {{1, 2}, {3, 6}, {7, 7}};
(*
  Location of constant-temperature boundaries
  among each list of corners
*)
domainXBathList["regular-long"] = {1};
domainXBathList["regular-short"] = {Way[patchedIntXList["regular-short"][[3]], 1]};
domainXBathList["irregular"] = {
  Way[patchedIntXList["irregular"][[1]], patchedIntXList["irregular"][[2]], 1/3],
  Way[patchedIntXList["irregular"][[2]], patchedIntXList["irregular"][[-2]], 1/5],
  Way[patchedCornerXList["irregular"][[-2]], patchedIntXList["irregular"][[-2]]]
};


(* ::Subsection:: *)
(*Fin example*)


(* ::Subsubsection:: *)
(*Physical constants*)


stefanBoltzmann = 5.67 * 10^-8 (* W /m^2 /K^4 *);


(*
  ## Emissivity ##
  Figure 4, in:
  Wade & Preedy (2003).
  Fujihokka: A high-emissivity approach to aluminum anodizing.
  Metal Finishing, 101(12), 8--13.
  <https://doi.org/10.1016/S0026-0576(03)90093-6>
*)
exampleEmissivity = 0.9;


(*
  ## Conductivity ##
  Figure 1, in:
  Cook, Moore, Matsumura, & Van der Meer. (1975)
  The Thermal and Electrical Conductivity of Aluminum.
  ORNL-5079.
  Oak Ridge National Lab.
  <https://doi.org/10.2172/5066461>
*)
exampleConductivity = 236 (* W /m /K *);


(*
  ## Radiation constant ##
  c == \[Epsilon] \[Sigma] / k
*)
With[
  {
    eps = exampleEmissivity,
    sigma = stefanBoltzmann,
    k = exampleConductivity,
    dummyForTrailingCommas = Null
  },
  exampleConstant = eps sigma / k (* /m /K^3 *);
];


(* ::Subsubsection:: *)
(*Interval*)


{exampleX1, exampleX2} = {0.5, 0.6};


(* ::Subsubsection:: *)
(*Profile*)


exampleYTraced[x_] := yTraLower[0][x] - yTraLower[0][exampleX1] // Evaluate;
exampleY2 = exampleYTraced[exampleX2];


(* ::Subsubsection:: *)
(*Base temperature T_2*)


With[
  {
    c = exampleConstant,
    x1Hat = exampleX1,
    x2Hat = exampleX2,
    dummyForTrailingCommas = Null
  },
  exampleT2[length_] :=
    x2Hat
    Divide[
      x2Hat - x1Hat,
      c * length
    ]^(1/3);
];


(* ::Subsubsection:: *)
(*Power per length p*)


With[
  {
    k = exampleConductivity,
    c = exampleConstant,
    x1Hat = exampleX1,
    x2Hat = exampleX2,
    y2Hat = exampleY2,
    dummyForTrailingCommas = Null
  },
  exampleP[length_] :=
    2 y2Hat k
    Divide[
      x2Hat - x1Hat,
      c * length
    ]^(1/3);
];


(* ::Subsubsection:: *)
(*Celsius, ugh*)


celsiusOffset = 273.15 (* K *);


(* ::Section:: *)
(*Traced boundary algebra*)


(* ::Subsection:: *)
(*General*)


With[{x = \[FormalX]},
  {
    {"dy/dx", yTraDer[x]},
    {"y", yTra[x]},
    {"y(1)", yTraMax},
    {"y(1)", yTraMax // N}
  } // TableForm
] // Ex["plane-traced-general.pdf"]


(* ::Subsection:: *)
(*Near x = 0*)


With[{x = \[FormalX]},
  Module[{yDerAsy, yAsy},
    yDerAsy = yTraDer[x] + O[x]^12;
    yAsy = Integrate[yDerAsy, x];
    {
      {"dy/dx", yDerAsy},
      {"y", yAsy}
    } // TableForm
  ]
] // Ex["plane-traced-x-near-0.pdf"]


(* ::Subsection:: *)
(*Near x = 1*)


With[{xi = \[FormalXi]},
  Module[{yDerAsy, yAsy},
    yDerAsy = -yTraDer[1 - xi] + O[xi]^(1/2);
    yAsy = Integrate[yDerAsy, xi];
    {
      {"dy/d\[Xi]", yDerAsy},
      {"y", yAsy}
    } // TableForm
  ]
] // Ex["plane-traced-x-near-1.pdf"]


(* ::Section:: *)
(*Traced boundary plot*)


(* ::Subsection:: *)
(*Without known solution*)


Plot[{yTra[x], -yTra[x]}, {x, 0, 1},
  AspectRatio -> Automatic,
  PlotRange -> Full,
  PlotStyle -> Black,
  PlotOptions[Axes] // Evaluate
] // Ex["plane-traced.pdf"]


(* ::Subsection:: *)
(*With known solution*)


Show[
  EmptyAxes[{0, 1}, {-yTraMax, yTraMax}],
  (* Known solution *)
  DensityPlot[x, {x, 0, 1}, {y, -yTraMax, yTraMax},
    ColorFunction -> "TemperatureMap",
    RegionFunction -> Function[{x, y}, Abs[y] < -yTra[x]]
  ],
  (* Traced boundaries *)
  Plot[{yTra[x], -yTra[x]}, {x, 0, 1},
    AspectRatio -> Automatic,
    PlotRange -> Full,
    PlotStyle -> Black,
    PlotOptions[Axes] // Evaluate
  ]
] // Ex["plane-traced-with-known.png"]


(* ::Section:: *)
(*Self-incident boundary ratio plots*)


(* ::Subsection:: *)
(*Exact expression*)


Module[
  {
    intervalList,
    numPoints,
    xSub, flatFractions, stringFractions,
    x1, x2,
    x1String, x2String, fileName,
    xValues, xRatioTable,
    dummyForTrailingCommas
  },
  (* List of intervals {x1, x2} to plot the ratio for *)
  intervalList = {
    {0, 1},
    {0, 1/2}, {1/2, 1},
    {0, 1/3}, {1/3, 2/3}, {2/3, 1},
    {0, 1/4}, {1/4, 1/2}, {1/2, 3/4}, {3/4, 1},
    Nothing
  };
  (* Number of points for sampling *)
  numPoints = 50;
  (* Nice labels *)
  xSub[n_] := Subscript[Italicise["x"], n];
  flatFractions = Rational[a_, b_] :> SeparatedRow["/"][a, b];
  stringFractions = {
    Rational[a_, b_] :> StringForm["``o``", a, b],
    n_Integer :> ToString[n]
  };
  (* Make plots *)
  Table[
    (* Get interval endpoints *)
    {x1, x2} = interval;
    (* Build table of values *)
    xValues = Subdivide[x1, x2, numPoints];
    xRatioTable = Table[{x, boundaryRatio[x1, x2][x]}, {x, xValues}];
    (* Plot and export *)
    {x1String, x2String} = {x1, x2} /. stringFractions;
    fileName = StringTemplate["plane-boundary-ratio-x1-``-x2-``.png"][x1String, x2String];
    ListPlot[xRatioTable
      , Joined -> True
      , AxesLabel -> Italicise /@ {"x", "R"}
      , PlotLabel -> {xSub[1], xSub[2]} == (interval /. flatFractions)
      , PlotOptions[Axes] // Evaluate
    ] // Ex[fileName]
    , {interval, intervalList}
  ] // Quiet
]


(* ::Subsection:: *)
(*Crude upper bound*)


Module[
  {
    intervalList,
    numPoints,
    xSub, flatFractions, stringFractions,
    x1, x2,
    x1String, x2String, fileName,
    xValues, xRatioTable,
    dummyForTrailingCommas
  },
  (* List of intervals {x1, x2} to plot the ratio for *)
  intervalList = {
    {0, 1},
    {0, 1/2}, {1/2, 1},
    {0, 1/3}, {1/3, 2/3}, {2/3, 1},
    {0, 1/4}, {1/4, 1/2}, {1/2, 3/4}, {3/4, 1},
    Nothing
  };
  (* Number of points for sampling *)
  numPoints = 50;
  (* Nice labels *)
  xSub[n_] := Subscript[Italicise["x"], n];
  flatFractions = Rational[a_, b_] :> SeparatedRow["/"][a, b];
  stringFractions = {
    Rational[a_, b_] :> StringForm["``o``", a, b],
    n_Integer :> ToString[n]
  };
  (* Make plots *)
  Table[
    (* Get interval endpoints *)
    {x1, x2} = interval;
    (* Build table of values *)
    xValues = Subdivide[x1, x2, numPoints];
    xRatioTable = Table[{x, boundaryRatio[x1, x2][x]}, {x, xValues}];
    (* Plot and export *)
    {x1String, x2String} = {x1, x2} /. stringFractions;
    fileName = StringTemplate["plane-boundary-ratio-with-bound-x1-``-x2-``.png"][x1String, x2String];
    Show[
      (* Crude upper bound *)
      Plot[boundaryRatioBound[x1, x2][x]
        , {x, x1, x2}
        , AxesLabel -> Italicise /@ {"x", "R"}
        , PlotLabel -> Column @ {
            {xSub[1], xSub[2]} == (interval /. flatFractions),
            Subscript[Italicise["R"], "ultra"] == boundaryRatioBoundUltra[x1, x2]
          }
        , PlotRange -> {{x1, x2}, {0, Automatic}}
        , PlotStyle -> Directive[Red, Dashed]
        , PlotOptions[Axes] // Evaluate
      ],
      (* Exact expression *)
      ListPlot[xRatioTable
        , Joined -> True
        , PlotRange -> All
      ],
      {}
    ] // Ex[fileName]
    , {interval, intervalList}
  ] // Quiet
]


(* ::Section:: *)
(*Fin example (length = 1 metre)*)


Module[
  {
    k, c, x1Hat, x2Hat, y2Hat,
    l12, w2, temp1, temp2, p,
    dummyForTrailingCommas
  },
  (* Abbreviations *)
  k = exampleConductivity;
  c = exampleConstant;
  {x1Hat, x2Hat} = {exampleX1, exampleX2};
  y2Hat = exampleY2;
  (* Computed quantities *)
  l12 = 1;
  w2 = 2 y2Hat l12 / (x2Hat - x1Hat);
  temp1 = x1Hat ((x2Hat - x1Hat) / (c l12))^(1/3);
  temp2 = x2Hat ((x2Hat - x1Hat) / (c l12))^(1/3);
  p = 2 y2Hat k ((x2Hat - x1Hat) / (c l12))^(1/3);
  (* Make table *)
  {
    {"Fin length", "L_12",
      l12 // Quantity[#, "Meters"] &
    },
    {"Fin width at base", "W_2",
      w2 // Quantity[#, "Meters"] &
    },
    {"Temperature at tip", "T_1",
      temp1 // Quantity[#, "Kelvins"] &,
      temp1 - celsiusOffset // Quantity[#, "DegreesCelsius"] &,
      {}
    },
    {"Temperature at base", "T_2",
      temp2 // Quantity[#, "Kelvins"] &,
      temp2 - celsiusOffset // Quantity[#, "DegreesCelsius"] &,
      temp2 == exampleT2[l12] (* Check *)
    },
    {"Power per length", "p",
      p // Quantity[#, "Watts" / "Meters"] &,
      {},
      p == exampleP[l12] (* Check *)
    },
    Nothing
  } // TableForm
] // Ex["plane-fin-1-metre.pdf"]


(* ::Section:: *)
(*Figure: Slab BVP (plane-slab-bvp.pdf)*)


Module[
  {
    xSlabMin, xSlabMax,
    ySlabMin, ySlabMax, ySlabMid,
    xRadiation, yRadiationMargin,
    numRadiation, yRadiationValues,
    radiationArrowStyle, radiationArrowSize,
    xAxisMin, xAxisMax, yAxis,
    axisStyle,
    textStyle,
    dummyForTrailingCommas
  },
  (* Slab dimensions *)
  {xSlabMin, xSlabMax} = {0, 1};
  {ySlabMin, ySlabMax} = {0, 1.7};
  ySlabMid = Way[ySlabMin, ySlabMax];
  (* Radiation arrow specs *)
  xRadiation = Way[xSlabMin, xSlabMax, -0.15];
  yRadiationMargin = Way[ySlabMin, ySlabMax, 0.05];
  numRadiation = 6;
  yRadiationValues = Subdivide[
    ySlabMin + yRadiationMargin,
    ySlabMax - yRadiationMargin,
    numRadiation - 1
  ];
  radiationArrowStyle = Arrowheads[0.025];
  radiationArrowSize = 2.2 Abs[xRadiation];
  (* x-axis specs *)
  xAxisMin = xRadiation - 2 radiationArrowSize;
  xAxisMax = xSlabMax + 2 radiationArrowSize;
  yAxis = ySlabMin - 0.7 radiationArrowSize;
  axisStyle = Directive[Arrowheads[0.05]];
  (* Make diagram *)
  textStyle = Style[#, LabelSize["Label"] - 1] & @* LaTeXStyle;
  Show[
    (* Slab horizontal boundaries *)
    Graphics @ {GeneralStyle["DefaultThick"],
      Line @ {{xSlabMin, ySlabMax}, {xSlabMax, ySlabMax}},
      Line @ {{xSlabMin, ySlabMin}, {xSlabMax, ySlabMin}},
      {}
    },
    (* Slab radiation boundary *)
    Graphics @ {BoundaryTracingStyle["Traced"],
      Line @ {{xSlabMin, ySlabMax}, {xSlabMin, ySlabMin}},
      {}
    },
    (* Radiation arrows *)
    Graphics @ {radiationArrowStyle,
      Table[
        SquigglyArrow[{xRadiation, yRadiation}, Pi, radiationArrowSize]
        , {yRadiation, yRadiationValues}
      ],
      Text[
        "radiation" // textStyle
        , {xRadiation - radiationArrowSize, ySlabMid}
        , {1.3, 0}
      ],
      {}
    },
    (* Constant temperature boundary *)
    Graphics @ {BoundaryTracingStyle["Contour"],
      Line @ {{xSlabMax, ySlabMax}, {xSlabMax, ySlabMin}},
      Text[
        Column[{"constant", "temperature"}
          , Alignment -> Center
          , Spacings -> 0
        ] // textStyle
        , {xSlabMax, ySlabMid}
        , {-1.2, 0}
      ],
      {}
    },
    (* x-axis *)
    Graphics @ {axisStyle,
      Arrow @ {{xAxisMin, yAxis}, {xAxisMax, yAxis}},
      Text[
        Italicise["x"] // textStyle
        , {xAxisMax, yAxis}
        , {-2.5, -0.2}
      ],
      {}
    },
    {}
    , ImageSize -> 0.45 ImageSizeTextWidth
  ]
] // Ex["plane-slab-bvp.pdf"]


(* ::Section:: *)
(*Figure: Traced boundaries, single spike (plane-traced-boundary-spike.pdf)*)


Module[{yMax, yMaxMore, xTerm},
  yMax = Ceiling[yTraMax, 0.1];
  yMaxMore = 1.2 yTraMax;
  xTerm = 1;
  Show[
    (* Pair of traced boundaries forming a spike *)
    Plot[{yTra[x], -yTra[x]}, {x, 0, 1},
      AspectRatio -> Automatic,
      ImageSize -> 240,
      LabelStyle -> LatinModernLabelStyle[12],
      PlotRange -> {-yMax, yMax},
      PlotStyle -> BoundaryTracingStyle["Traced"],
      PlotOptions[Axes] // Evaluate
    ],
    (* Critical terminal curve *)
    Graphics @ {BoundaryTracingStyle["Contour"],
      Line @ {{xTerm, -yMaxMore}, {xTerm, yMaxMore}}
    }
  ]
] // Ex["plane-traced-boundary-spike.pdf"]


(* ::Section:: *)
(*Figure: Traced boundaries, various (plane-traced-boundary-various.pdf)*)


(*
  REMEMBER:
    the positive signed branch is the lower;
    the negative signed branch is the upper.
  This is because the sign in the tracing equation is \[MinusPlus].
  Therefore in a spike (shaped like <),
  the lower-branch curve is actually higher.
  However, writing each curve as y == c \[MinusPlus] 2F1(x; ...),
  note that the lower-branch curve has lower c.
  Now, to generate assorted patchings of traced boundaries,
  we choose the corners (x_1, y_1), ..., (x_n, y_n)
  with y_1 < ... < y_n
  and determine the intersections
  between lower-branch(i) and upper-branch(i+1)
  for i == 1, ..., n - 1.
  The corners (x_i, y_i) cannot be chosen arbitrarily,
  so don't be dumb.
 *)
Module[
 {yTraUpper, yTraLower,
  cornerListList,
  xMin, xMax, yMin, yMax,
  imageSize,
  plotList,
  xCornerList, yCornerList,
  cUpperList, cLowerList,
  xCorner, yCorner,
  cUpper, cLower,
  n, xIntList, yIntList,
  restrict,
  xInt, yInt,
  xLeft, xRight,
  xTerm
 },
  (* Avoid confusion (note yTra == -2F1(x; ...)) *)
  yTraUpper[c_][x_] := c + yTra[x];
  yTraLower[c_][x_] := c - yTra[x];
  (* List of lists of corners *)
  cornerListList = List[
    Table[{0.4, y}, {y, Subdivide[-0.6, 0.6, 4]}],
    Table[{0.8, y}, {y, Subdivide[-0.7, 0.7, 8]}],
    {0.1, 0.3} # & /@ {
      {4, -2.2},
      {9, -1.3},
      {7.7, -0.5},
      {2.3, -0.2},
      {4.5, 0},
      {7, 0.5},
      {3, 2}
    }
  ];
  (* Critical terminal curve *)
  xTerm = 1;
  (* Plot range *)
  xMin = 0.2;
  xMax = 1.05 xTerm;
  yMax = 0.7;
  imageSize = 240;
  (* Build a plot for each of these lists *)
  plotList = Table[
    {xCornerList, yCornerList} = Transpose[cornerList];
    (*
      List of constants {cUpper, cLower}
      for the pair of traced boundaries
      through each corner point
     *)
    {cUpperList, cLowerList} =
      Transpose @ Table[
        {xCorner, yCorner} = corner;
        cUpper = yCorner - yTra[xCorner];
        cLower = yCorner + yTra[xCorner];
        {cUpper, cLower}
      , {corner, cornerList}];
    (*
      List of intersections between
      lower-branch(i) and upper-branch(i+1)
     *)
    n = Length[cornerList];
    {xIntList, yIntList} =
      Transpose @ Table[
        cLower = cLowerList[[i]];
        cUpper = cUpperList[[i + 1]];
        xInt = SeekRoot[
          yTraUpper[cUpper][#] - yTraLower[cLower][#] &,
          {0, 1}
        ];
        yInt = yTraUpper[cUpper][xInt];
        {xInt, yInt}
      , {i, n - 1}];
    (* Restrict domain *)
    restrict[x0_, x1_] :=
      Piecewise[
        {{1, x0 <= # <= x1}},
        Indeterminate
      ] &;
    (* Plot *)
    Show[
      (* Traced boundaries *)
      Plot[
        Table[
          {
            (* Upper-branch(i) *)
            xLeft = xCornerList[[i]];
            xRight = If[i > 1,
              xIntList[[i - 1]],
              Max[xIntList]
            ];
            cUpper = cUpperList[[i]];
            (
              yTraUpper[cUpper][x] restrict[xLeft, xRight][x]
            ),
            (* Lower-branch(i) *)
            xLeft = xCornerList[[i]];
            xRight = If[i < n,
              xIntList[[i]],
              Max[xIntList]
            ];
            cLower = cLowerList[[i]];
            (
              yTraLower[cLower][x] restrict[xLeft, xRight][x]
            )
          }
        , {i, n}] // Evaluate,
        {x, 0, 1},
        AspectRatio -> Automatic,
        Axes -> None,
        ImageSize -> 240,
        PlotRange -> {{xMin, xMax}, {-yMax, yMax}},
        PlotStyle -> BoundaryTracingStyle["Traced"]
      ],
      (* Critical terminal curve *)
      Graphics @ {BoundaryTracingStyle["Contour"],
        Line @ {{xTerm, -yMax}, {xTerm, yMax}}
      }
    ]
  , {cornerList, cornerListList}]
  // GraphicsRow[#,
    Spacings -> {
      {0, -0.2, 0.55} imageSize,
      Automatic
    }
  ] &
] // Ex["plane-traced-boundary-various.pdf"]


(* ::Section:: *)
(*Figure: Traced boundaries (plane-traced-boundaries.pdf)*)


Module[
 {xMin, xMax, yMax,
  xTerm,
  cMax, cList
 },
  (* Plot range *)
  xMin = 0;
  xMax = 1;
  yMax = Ceiling[2 yTraMax, 0.1];
  (* Critical terminal curve *)
  xTerm = 1;
  (* Values of integration constant *)
  cMax = 0.95 yMax;
  cList = Subdivide[-yMax, yMax, 8];
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      AspectRatio -> Automatic,
      FrameLabel -> {
        Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
        Italicise["y"]
      },
      FrameTicksStyle -> LabelSize["Tick"],
      ImageSize -> 0.45 ImageSizeTextWidth,
      LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    ],
    (* Traced boundaries *)
    Table[
      Plot[c + yTra[x] {1, -1}, {x, xMin, xMax},
        PlotPoints -> 3,
        PlotRange -> {-yMax, yMax},
        PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
    , {c, cList}],
    (* Critical terminal curve *)
    ParametricPlot[
      {xTerm, y}, {y, -yMax, yMax},
      PlotPoints -> 2,
      PlotStyle -> BoundaryTracingStyle["Traced"]
    ]
  ]
] // Ex["plane-traced-boundaries.pdf"]


(* ::Section:: *)
(*Figure: Traced boundaries, patched (plane-traced-boundaries-patched.pdf)*)


Module[
 {xTerm,
  xMin, xMax, yMin, yMax,
  plotList,
  plotPointsGeneral, plotPointsPatched,
  textStyle,
  n, cUpperList, cLowerList, xCornerList, xIntList,
  cUpper, cLower,
  xLeft, xRight
 },
  (* Critical terminal curve *)
  xTerm = 1;
  (* Plot range *)
  xMin = 0.2;
  xMax = 1.1 xTerm;
  yMax = 0.7;
  (* Plot points *)
  plotPointsGeneral = 3;
  plotPointsPatched = 2;
  (* Styles *)
  textStyle = Style[#, LabelSize["Straight"]] & @* LaTeXStyle;
  (* Build a plot for each list of corners *)
  plotList = Table[
    (* *)
    n = patchedCornerNum[id];
    cUpperList = patchedCUpperList[id];
    cLowerList = patchedCLowerList[id];
    xCornerList = patchedCornerXList[id];
    xIntList = patchedIntXList[id];
    (* Plot *)
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        AspectRatio -> Automatic,
        Frame -> None,
        ImageSize -> Automatic
      ],
      (* General boundaries: upper-branch(i) *)
      Table[
        cUpper = cUpperList[[i]];
        Plot[
          yTraUpper[cUpper][x],
          {x, 0, 1},
          PlotPoints -> plotPointsGeneral,
          PlotRange -> {-yMax, yMax},
          PlotStyle -> BoundaryTracingStyle["Background"]
        ]
      , {i, n}],
      (* General boundaries: lower-branch(i) *)
      Table[
        cLower = cLowerList[[i]];
        Plot[
          yTraLower[cLower][x],
          {x, 0, 1},
          PlotPoints -> plotPointsGeneral,
          PlotRange -> {-yMax, yMax},
          PlotStyle -> BoundaryTracingStyle["Background"]
        ]
      , {i, n}],
      (* Patched portions: upper-branch(i) *)
      Table[
        xLeft = xCornerList[[i]];
        xRight = If[i > 1,
          xIntList[[i - 1]],
          Max[xIntList]
        ];
        cUpper = cUpperList[[i]];
        Plot[
          yTraUpper[cUpper][x],
          {x, xLeft, xRight},
          PlotPoints -> plotPointsPatched,
          PlotRange -> {-yMax, yMax},
          PlotStyle -> BoundaryTracingStyle["Traced"]
        ]
      , {i, n}],
      (* Patched portions: lower-branch(i) *)
      Table[
        xLeft = xCornerList[[i]];
        xRight = If[i < n,
          xIntList[[i]],
          Max[xIntList]
        ];
        cLower = cLowerList[[i]];
        Plot[
          yTraLower[cLower][x],
          {x, xLeft, xRight},
          PlotPoints -> plotPointsPatched,
          PlotRange -> {-yMax, yMax},
          PlotStyle -> BoundaryTracingStyle["Traced"]
        ]
      , {i, n}],
      (* Critical terminal curve *)
      ParametricPlot[
        {xTerm, y}, {y, -yMax, yMax},
        PlotPoints -> 2,
        PlotStyle -> BoundaryTracingStyle["Terminal"]
      ],
      Graphics @ {
        Text[
          Italicise["x"] == 1
          , {1, 0}
          , {0, 1}
          , {0, 1}
        ] // textStyle
      },
      {}
    ]
  , {id, patchedIdList}]
  // GraphicsRow[#,
    Spacings -> {
      0.3 ImageSizeTextWidth,
      Automatic
    }
  ] &
  // Show[#, ImageSize -> ImageSizeTextWidth] &
] // Ex["plane-traced-boundaries-patched.pdf"]


(* ::Section:: *)
(*Figure: Domains (plane-domains.pdf)*)


(* ::Subsection:: *)
(*Domains*)


Module[
 {xTerm,
  xMin, xMax, yMin, yMax,
  plotPointsPatched,
  textStyle, textStyleLabel, textStyleBracket,
  plotList,
  n, cUpperList, cLowerList, xCornerList, xIntList,
  iRangeList, xBathList,
  iMin, iMax, xBath, yBathBottom, yBathTop,
  cUpper, cLower,
  xLeft, xRight,
  domainLabel,
  dummyForTrailingCommas
 },
  (* Critical terminal curve *)
  xTerm = 1;
  (* Plot range *)
  xMin = 0.2;
  xMax = 1.1 xTerm;
  yMax = 0.7;
  (* Plot points *)
  plotPointsPatched = 2;
  (* Styles *)
  textStyle = Style[#, LabelSize["Straight"]] & @* LaTeXStyle;
  textStyleLabel = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  (* Build a plot for each list of corners *)
  plotList = Table[
    (* *)
    n = patchedCornerNum[id];
    cUpperList = patchedCUpperList[id];
    cLowerList = patchedCLowerList[id];
    xCornerList = patchedCornerXList[id];
    xIntList = patchedIntXList[id];
    iRangeList = domainCornerRangeList[id];
    xBathList = domainXBathList[id];
    (* Plot *)
    Show[
      EmptyAxes[{xMin, xMax}, {-yMax, yMax},
        AspectRatio -> Automatic,
        Axes -> None,
        ImageSize -> Automatic
      ],
      (* Critical terminal curve *)
      ParametricPlot[
        {xTerm, y}, {y, -yMax, yMax},
        PlotPoints -> 2,
        PlotStyle -> BoundaryTracingStyle["Terminal", "Background"]
      ],
      Graphics @ {Gray,
        Text[
          Italicise["x"] == 1
          , {1, 0}
          , {0, 1}
          , {0, 1}
        ] // textStyle
      },
      (* Domains *)
      Table[
        {iMin, iMax} = iRangeList[[j]];
        {
          (* Constant-temperature (Dirichlet) boundary *)
          xBath = xBathList[[j]];
          yBathBottom = yTraUpper[cUpperList[[iMin]]] @ xBath;
          yBathTop = yTraLower[cLowerList[[iMax]]] @ xBath;
          ParametricPlot[
            {xBath, y}, {y, yBathBottom, yBathTop},
            PlotPoints -> 2,
            PlotStyle -> BoundaryTracingStyle["Contour"]
          ],
          (* Patched portions: upper-branch(i) *)
          Table[
            xLeft = xCornerList[[i]];
            xRight = If[i > iMin, xIntList[[i - 1]], xBath];
            cUpper = cUpperList[[i]];
            Plot[
              yTraUpper[cUpper][x],
              {x, xLeft, xRight},
              PlotPoints -> plotPointsPatched,
              PlotRange -> {-yMax, yMax},
              PlotStyle -> BoundaryTracingStyle["Traced"]
            ]
          , {i, iMin, iMax}],
          (* Patched portions: lower-branch(i) *)
          Table[
            xLeft = xCornerList[[i]];
            xRight = If[i < iMax, xIntList[[i]], xBath];
            cLower = cLowerList[[i]];
            Plot[
              yTraLower[cLower][x],
              {x, xLeft, xRight},
              PlotPoints -> plotPointsPatched,
              PlotRange -> {-yMax, yMax},
              PlotStyle -> BoundaryTracingStyle["Traced"]
            ]
          , {i, iMin, iMax}]
        }
      , {j, Length[iRangeList]}],
      (* Domain labels *)
      domainLabel[str_][scaledX_, scaledY_] :=
        Text[
          Row @ {
            "(" // textStyleBracket,
            str,
            ")" // textStyleBracket,
          } // textStyleLabel
          , {Way[xMin, xMax, scaledX], scaledY * yMax}
        ];
      Graphics @ {
        Which[
          id == patchedIdList[[1]],
            domainLabel["a"][0.4, 0],
          id == patchedIdList[[2]],
            domainLabel["b"][0.5, 0],
          id == patchedIdList[[3]],
            {
              domainLabel["c"][0.3, 0.75],
              domainLabel["d"][0.11, 0.05],
              domainLabel["e"][0.5, -0.62],
              Nothing
            },
          True, {}
        ]
      },
      {}
    ]
  , {id, patchedIdList}];
  (* Combine *)
  GraphicsRow[plotList
    , ImageSize -> ImageSizeTextWidth
    , Spacings -> {0.3 ImageSizeTextWidth, Automatic}
  ]
] // Ex["plane-domains.pdf"]


(* ::Subsection:: *)
(*Legend*)


Module[{legendCurves},
  legendCurves =
    CurveLegend[
      BoundaryTracingStyle @* ReleaseHold /@
        {"Traced", "Contour", Hold @ Sequence["Terminal", "Background"]},
      {"radiation", "constant temperature", "critical terminal curve"}
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Legend"]
    ];
  GraphicsGrid[{legendCurves}
    , Alignment -> Left
    , ImageSize -> ImageSizeTextWidth
    , ItemAspectRatio -> 0.11
    , Spacings -> {{0, -0.3, 0.1} ImageSizeTextWidth, 0}
  ]
] // Ex["plane-domains-legend.pdf"]


(* ::Section:: *)
(*Figure: Self-viewing radiation ratio (plane-self-viewing-ratio.pdf)*)


Module[
  {
    labelPlacementFromInterval, intervalList,
    numPoints, data,
    x1, x2, xValues,
    textStyle, textStyleBracket,
    intervalText,
    proportion, offset, xLabel, rLabel,
    dummyForTrailingCommas
  },
  (* List of intervals {x1, x2} to plot the ratio for *)
  labelPlacementFromInterval = Association[
    {0, 0.3} -> {1, {-1.2, 0}},
    {0.15, 0.45} -> {1, {-0.9, 0.8}},
    {0.3, 0.6} -> {1, {-0.9, 0.8}},
    {0.45, 0.75} -> {1, {-1.05, 0.4}},
    {0.6, 0.9} -> {0, {0.6, -0.8}},
    {0.7, 1} -> {0.3, {0, -1.4}},
    Nothing
  ];
  intervalList = Keys[labelPlacementFromInterval];
  (* Number of points for sampling *)
  numPoints = 50;
  (* Compute data points to plot *)
  Table[
    {x1, x2} = interval;
    xValues = Subdivide[x1, x2, numPoints];
    If[x1 == 0,
      xValues = Join[
        Subdivide[x1, xValues[[4]], 20],
        xValues[[5 ;;]]
      ];
    ];
    If[x2 == 1,
      xValues = Join[
        xValues[[;; -5]],
        Subdivide[xValues[[-4]], x2, 20]
      ];
    ];
    data[interval] =
      Cases[
        Table[{x, boundaryRatio[x1, x2][x]}, {x, xValues}],
        {_?NumericQ, _?NumericQ}
      ] // Quiet;
    , {interval, intervalList}
  ];
  (* Text style *)
  textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  intervalText[interval_] :=
    Text[
      {x1, x2} = interval;
      {proportion, offset} = labelPlacementFromInterval[interval];
      Row @ {
        "(" // textStyleBracket,
        "\[NegativeVeryThinSpace]",
        x1,
        ",\[ThinSpace]",
        x2,
        ")" // textStyleBracket
      } // textStyle
      ,
        xLabel = Way[x1, x2, proportion];
        rLabel = Interpolation[data[interval], xLabel];
        {xLabel, Log[rLabel]}
      , offset
    ];
  (* Make plot *)
  ListLogPlot[
    Table[data[interval], {interval, intervalList}]
    , AxesLabel -> {
        Italicise["x"] // Margined @ {{0, 1}, {5, 0}},
        Italicise["R"]
      }
    , Epilog -> {
        {
          BoundaryTracingStyle["Contour"],
          Line @ {{0, #}, {1, #}} & [1/100 // Log]
        },
        Table[intervalText[interval], {interval, intervalList}],
        {}
      }
    , ImageSize -> 0.55 ImageSizeTextWidth
    , Joined -> True
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotRange -> {{0, 1}, {All, 10}}
    , PlotRangeClipping -> False
    , PlotStyle -> Black
    , TicksStyle -> LabelSize["Tick"]
  ]
] // Ex["plane-self-viewing-ratio.pdf"]


(* ::Section:: *)
(*Figure: example fin base temperature (plane-fin-base-temperature)*)


Plot[
  exampleT2[l12] - celsiusOffset
  , {l12, 0, 6}
  , AxesLabel -> {
      (* L_12 / m *)
      SeparatedRow["VeryThin"][
        Subscript[Italicise["L"], "12"],
        Style["/", Magnification -> 1.25],
        "m"
      ] // Margined @ {{-2, 0}, {2, 0}},
      (* t_2 / \[Degree]C *)
      SeparatedRow["VeryThin"][
        Subscript[Italicise["t"], "2"],
        Style["/", Magnification -> 1.25],
        SeparatedRow[][Style["\[Degree]", Magnification -> 1.25], "C"]
      ] // Margined @ {{0, 0}, {-3, -5}}
    }
  , ImageSize -> 0.5 ImageSizeTextWidth
  , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
  , PlotStyle -> Black
  , TicksStyle -> LabelSize["Tick"]
] // Ex["plane-fin-base-temperature.pdf"]
