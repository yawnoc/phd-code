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


yTraDer = Function[{x}, -x^4 / Sqrt[1 - x^8]];
yTra = Function[{x}, Integrate[yTraDer[x], x] // Evaluate];
yTraMax = Abs @ yTra[1];


With[{x = \[FormalX]},
  {yTraDer[x], yTra[x], yTraMax}
]


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
