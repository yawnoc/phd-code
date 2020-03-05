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


(* ::Subsection:: *)
(*Known solution*)


(* ::Subsubsection:: *)
(*T*)


(* ::Text:: *)
(*See (r5.10) (Page r5-2).*)


tKnown[b_][x_, y_] := 1 - b Cos[x] Cosh[y];


(* ::Subsubsection:: *)
(*Location of straight contour*)


xStraight = Pi/2;


tKnown[\[FormalCapitalB]][xStraight, \[FormalY]] == 1


(* ::Subsubsection:: *)
(*P = \[PartialD]T/\[PartialD]x, Q = \[PartialD]T/\[PartialD]y*)


(* ::Text:: *)
(*See (r5.13) and (r5.14) (Page r5-2).*)


p[b_][x_, y_] := D[tKnown[b][x, y], x] // Evaluate;
q[b_][x_, y_] := D[tKnown[b][x, y], y] // Evaluate;


With[{x = \[FormalX], y = \[FormalY], b = \[FormalCapitalB]},
  {
    p[b][x, y] == b Sin[x] Cosh[y],
    q[b][x, y] == -b Cos[x] Sinh[y]
  }
]


(* ::Subsection:: *)
(*Flux function F*)


(* ::Text:: *)
(*See (r5.15) (Page r5-2).*)


f[a_, b_][x_, y_] := -tKnown[b][x, y]^4 / a // Evaluate;


(* ::Subsection:: *)
(*Viability \[CapitalPhi]*)


(* ::Text:: *)
(*See (r5.17) (Page r5-2).*)


vi[a_, b_][x_, y_] := p[b][x, y]^2 + q[b][x, y]^2 - f[a, b][x, y]^2 // Evaluate;


(* ::Subsection:: *)
(*Simple case (B = 1)*)


(* ::Subsubsection:: *)
(*Critical terminal points along y = 0*)


(* ::Text:: *)
(*Observe that \[CapitalPhi](y = 0) = (1 - c^2) - (1 - c)^8 / A^2, where c = cos(x).*)
(*The critical terminal point at x = 0 is trivial.*)


x0Simp[a_] :=
  With[{c = \[FormalC]},
    ArcCos[c]
      /. Last @ Solve[(1 - c^2) - (1 - c)^8 / a^2 == 0, c, Reals]
  ] // Evaluate;


(* ::Subsubsection:: *)
(*Representative values of A*)


aValuesSimp = {1/10, 1, 2};


aNamesSimp = AssociationThread[
  aValuesSimp,
  {"small", "med", "large"}
];


aValuesSimpConvex = {8/10, 9/10};


aNamesSimpConvex = AssociationThread[
  aValuesSimpConvex,
  {"convex_8", "convex_9"}
];


(* ::Subsection:: *)
(*Starting points for boundary tracing*)


(* ::Subsubsection:: *)
(*System of ODES for terminal curve \[CapitalPhi] = 0*)


viContourSystem[a_, b_, sSign_: 1] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{fun, p, q, slope},
      (* Scalar function \[CapitalPhi] whose contours are sought *)
      fun = vi[a, b][x, y];
      (* Components of the gradient vector *)
      p = D[fun, x];
      q = D[fun, y];
      (* Magnitude of the gradient vector *)
      slope = Sqrt[p^2 + q^2];
      (* Return system of ODEs *)
      {
        x' == q / slope,
        y' == -p / slope
      } /. {
        x' -> Sign[sSign] x'[s],
        y' -> Sign[sSign] y'[s],
        x -> x[s],
        y -> y[s],
        List -> Sequence
      }
    ]
  ] // Evaluate;


(* ::Subsubsection:: *)
(*System of ODES for T-contours*)


xyContourSystem[b_] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    With[
     {p = p[b][x, y],
      q = q[b][x, y]
     },
      Module[{slope, xDer, yDer},
        (* Magnitude of the gradient vector *)
        slope = Sqrt[p^2 + q^2];
        (* Return system of ODEs *)
        {
          x' == q / slope,
          y' == -p / slope
        } /. {
          x' -> x'[s],
          y' -> y'[s],
          x -> x[s],
          y -> y[s]
        }
      ]
    ]
  ] // Evaluate;


(* ::Subsubsection:: *)
(*Simple case (B = 1)*)


(* Non-trivial critical terminal point (x_0, 0), which is hyperbolic *)
Table[
  startXYSimp[a]["hyperbolic"] = {
    {x0Simp[a], 0}
  };
, {a, aValuesSimp}];


(* Starting points along contour through (1/2 x_0, 0) *)
Table[
  startXYSimp[a]["contour"] =
    With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
      Module[{b, nMax, yMax, sMax, xyContour},
        b = 1;
        nMax = 8;
        yMax = 2;
        (* (Probable) upper bound for arc length traversed *)
        sMax = 3/2 * yMax;
        (* Terminal curve *)
        xyContour =
          NDSolveValue[
            {
              xyContourSystem[b],
              x[0] == 1/2 x0Simp[a], y[0] == 0,
              WhenEvent[Abs @ y[s] > yMax, "StopIntegration"]
            }, {x, y}, {s, -sMax, sMax},
            NoExtrapolation
          ];
        (* Actual arc length traversed *)
        sMax = DomainEnd[xyContour];
        (* Starting points along the terminal curve *)
        Table[
          xyContour[s] // Through // Rationalize[#, 0] &
        , {s, Subdivide[-sMax, sMax, nMax]}]
      ]
    ];
, {a, aValuesSimp}];


(* ::Subsection:: *)
(*Terminal curve x = x(s), y = y(s)*)


(* ::Subsubsection:: *)
(*Simple case (B = 1)*)


xyTermSimp[a_?NumericQ] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{b, sMax},
      b = 1;
      sMax = 5;
      NDSolveValue[
        {
          viContourSystem[a, b],
          x[0] == x0Simp[a], y[0] == 0
        }, {x, y}, {s, -sMax, sMax},
        NoExtrapolation
      ]
    ]
  ];


(* ::Subsection:: *)
(*Traced boundaries x = x(y)*)


(* ::Subsubsection:: *)
(*Traced boundary derivative dx/dy*)


(* ::Text:: *)
(*See (r5.18) (Page r5-3).*)
(*The lower branch is used since this corresponds*)
(*to the candidate boundary for y > 0.*)


xTraDer[a_, b_] := Function[{x, y},
  With[
   {p = p[b][x, y],
    q = q[b][x, y],
    f = f[a, b][x, y],
    vi = vi[a, b][x, y]
   },
    Divide[
      p q + f Sqrt[vi],
      q^2 - f^2
    ]
  ] // Evaluate
] // Evaluate;


(* ::Subsubsection:: *)
(*Traced boundary curvature*)


(* ::Text:: *)
(*See (r5.20) (Page r5-5).*)


curTra[a_, b_] := Function[{x, y},
  Module[{d, xDer, xDer2, cur},
    (* Abbreviation for y-derivative *)
    d = Dt[#, y, Constants -> {a, b}] &;
    (* x' and x'' *)
    xDer = xTraDer[a, b][x, y];
    xDer2 = d[xDer];
    (* Curvature evaluated *)
    cur = xDer2 /. {d[x] -> xDer}
  ] // Evaluate
] // Evaluate;


(* ::Subsubsection:: *)
(*Simple case (B = 1)*)


(* ::Subsubsubsection:: *)
(*Candidate boundary*)


xTraCandSimp[a_?NumericQ, terminateAtStraightContour_: False] :=
  With[{x = \[FormalX]},
    Module[{b, yMax},
      b = 1;
      yMax = 5;
      NDSolveValue[
        {
          x'[y] == Re @ xTraDer[a, b][x[y], y],
          x[0] == x0Simp[a],
          WhenEvent[
            terminateAtStraightContour && x[y] > xStraight,
            "StopIntegration"
          ]
        }, x, {y, 0, yMax},
        NoExtrapolation
      ]
    ]
  ];


(* ::Subsubsubsection:: *)
(*Critical value y = y(A) for inflection along x = \[Pi]/2.*)


(* ::Text:: *)
(*See (r5.24) (Page r5-6).*)


yCurCritSimp[a_] :=
  With[{y = \[FormalY], s = \[FormalCapitalS]},
    ArcSinh[s] /. First @ Solve[
      {
        2 s - (1 + s^2) (a^2 s + 4 Sqrt[a^2 (1 + s^2) - 1]) == 0,
        0 < a < 1
      },
      s, Reals
    ]
  ] // Evaluate;


yCurCritSimp[\[FormalCapitalA]]


(* ::Subsubsubsection:: *)
(*A_i (inflection dimensionless group)*)


(* ::Text:: *)
(*See (r5.25) (Page r5-6).*)


(* Compute A_i using the bisection algorithm *)
(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
aInflSimp = Module[{dest, aMin, aMax, a, num},
  dest = "cosine_simple-a-inflection.txt";
  ExportIfNotExists[dest,
    (* A_i is around 0.8 *)
    aMin = 7/10;
    aMax = 9/10;
    (* Solve x(y(A)) - \[Pi]/2 == 0 *)
    SeekRootBisection[
      xTraCandSimp[#] @ yCurCritSimp[#] - Pi / 2 &,
      {aMin, aMax},
      "ReturnIterations" -> True
    ] // Compress
  ];
  (* Import *)
  {a, num} = Import[dest] // Uncompress;
  (* Print iterations used *)
  Print @ FString["Bisection algorithm: {num} iterations"];
  (* Return A_i *)
  a
]


(* ::Subsubsubsection:: *)
(*Approximate determination of A_i*)


(* ::Text:: *)
(*See (r5.29) (Page r5-7).*)


aInflSimpApprox =
  With[{a = \[FormalA]},
    a /. First @ Solve[
      {
        a^6 - a^4 + 44 a^2 - 28 == 0,
        a > 0
      },
      a, Reals
    ]
  ];
aInflSimpApprox // N


(* ::Subsubsubsection:: *)
(*Convex domain aspect ratio*)


(* ::Text:: *)
(*Ratio of height to width of the thin, lens-like convex domains.*)


aspectRatioSimp[a_] :=
  Module[{xRad, xStart, xEnd, yStart, yEnd},
    (* Traced boundary *)
    xRad = xTraCandSimp[a, True];
    xStart = xRad[0];
    xEnd = xStraight;
    yStart = 0;
    yEnd = DomainEnd[xRad];
    (* Return aspect ratio *)
    2 (yEnd - yStart) / (xEnd - xStart)
  ];


(* ::Subsection:: *)
(*Traced boundaries x = x(s), y = y(s)*)


(* ::Subsubsection:: *)
(*System of ODES for tracing*)


xyTraSystem[a_, b_] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    With[
     {p = p[b][x, y],
      q = q[b][x, y],
      f = f[a, b][x, y],
      vi = vi[a, b][x, y]
     },
      Module[{grad2, xDer, yDer},
        (* Square of gradient *)
        grad2 = p^2 + q^2;
        (* Return system of ODEs *)
        xDer = (-q f + p Re @ Sqrt[vi]) / grad2;
        yDer = (+p f + q Re @ Sqrt[vi]) / grad2;
        {x' == xDer, y' == yDer} /. {
          x' -> x'[s],
          y' -> y'[s],
          x -> x[s],
          y -> y[s]
        }
      ]
    ]
  ] // Evaluate;


(* ::Subsubsection:: *)
(*Simple case (B = 1)*)


(* ::Subsubsubsection:: *)
(*Generic boundaries (for the chosen representative values of A)*)


Module[
 {b,
  idList,
  sMax,
  xyInitList, viTol,
  xInit, yInit
 },
  b = 1;
  idList = {"contour", "hyperbolic"};
  sMax = 4;
  Table[
    Table[
      xyInitList = startXYSimp[a][id];
      viTol = 10^-6;
      (*
        NOTE: leniency is required;
        when travelling away from y == 0,
        the traced boundaries asymptote very quickly
        towards the terminal curve.
       *)
      xyTraSimp[a][id] =
        Table[
          {xInit, yInit} = xyInit;
          With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
            NDSolveValue[
              {
                xyTraSystem[a, b],
                x[0] == xInit, y[0] == yInit,
                WhenEvent[
                  {
                    vi[a, b][x[s], y[s]] < viTol,
                    tKnown[b][x[s], y[s]] < 0
                  },
                  "StopIntegration"
                ]
              }, {x, y}, {s, -sMax, sMax}
            ]
          ]
        , {xyInit, xyInitList}];
    , {id, idList}]
  , {a, aValuesSimp}]
];


(* ::Subsubsubsection:: *)
(*Candidate boundary*)


(*
  "Cand" means candidate,
  referring to the traced boundary through (x_0, 0)
  which is almost the same as the terminal curve.
  The upper branch is taken for y > 0;
  the lower branch is taken for y < 0.
 *)
xyTraCandSimp[a_?NumericQ] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{b, sMax},
      b = 1;
      sMax = 5;
      NDSolveValue[
        {
          xyTraSystem[a, b],
          x[0] == x0Simp[a], y[0] == 0
        }, {x, y}, {s, 0, sMax},
        NoExtrapolation
      ]
    ]
  ];


(* ::Subsection:: *)
(*Italicised symbols*)


aIt = Italicised["A"];
bIt = Italicised["B"];
sIt = Italicised["s"];
xIt = Italicised["x"];
yIt = Italicised["y"];


(* ::Subsection:: *)
(*Options for exported GIFS*)


gifOpts = Sequence[
  AnimationRepetitions -> Infinity,
  "DisplayDurations" -> 0.5
];


(* ::Subsection:: *)
(*Global styles for plots*)


textStyle = Style[#, 18] &;


contStyle = LightGray;
straightStyle = Directive[Dashed, Magenta];
nonStyle = Directive[Opacity[0.7], LightGray];
termStyle = Pink;
unphysStyle = Black;


upperStyle = Blue;
lowerStyle = Red;
convexStyle = Black;


simp0Style = Blue;


pointStyle = PointSize[Large];
glowStyle = Directive[Thick, Yellow, Opacity[0.7]];


inflDotStyle = Directive[Red, Opacity[0.7], pointStyle];


(* ::Section:: *)
(*Known solution*)


(* ::Subsection:: *)
(*Algebra check for \[CapitalPhi]*)


(* ::Text:: *)
(*See (r5.17) (Page r5-2).*)


With[{a = \[FormalCapitalA], b = \[FormalCapitalB], x = \[FormalX], y = \[FormalY]},
  vi[a, b][x, y] ==
    b^2 (Sin[x]^2 + Sinh[y]^2)
    - (1 - b Cos[x] Cosh[y])^8 / a^2
    // FullSimplify
]


(* ::Subsection:: *)
(*Interactive visualiser for known solution*)


DynamicModule[
 {bInit, bMin, bMax,
  xMin, xMax, yMax
 },
  (* Values of B *)
  bInit = 1;
  bMin = 0;
  bMax = 5;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2;
  yMax = 3;
  (* Plot *)
  Manipulate[
    Plot3D[
      tKnown[b][x, y],
      {x, xMin, xMax}, {y, -yMax, yMax},
      ClippingStyle -> unphysStyle,
      PlotLabel -> BoxedLabel[bIt == N[b]],
      PlotOptions[Axes] // Evaluate,
      PlotRange -> {0, Automatic}
    ]
  , {{b, bInit, bIt}, bMin, bMax}]
]


(* ::Subsection:: *)
(*Interactive visualiser for unphysical domain (T < 0)*)


DynamicModule[
 {bInit, bMin, bMax,
  xMin, xMax, yMax,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinCont, xMaxCont, yMaxCont,
  numTo1, numBeyond1
 },
  (* Values of B *)
  bInit = 1;
  bMin = 0;
  bMax = 5;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 5/4;
  yMax = 2;
  (* Margin *)
  eps = 0.05;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin - eps;
  xMaxUnphys = xMax + eps;
  yMaxUnphys = yMax + eps;
  (* Plot range for contours *)
  xMinCont = xMin - eps;
  xMaxCont = xMax + eps;
  yMaxCont = yMax + eps;
  (* Number of contours *)
  numTo1 = 5;
  numBeyond1 = 3;
  (* Plot *)
  Manipulate[
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        ImageSize -> 240,
        PlotLabel -> BoxedLabel[bIt == N[b]]
      ],
      (* Unphysical domain *)
      RegionPlot[
        tKnown[b][x, y] < 0,
        {x, xMin, xMax}, {y, -yMax, yMax},
        BoundaryStyle -> unphysStyle,
        PlotPoints -> 50,
        PlotStyle -> unphysStyle
      ],
      (* Known solution contours *)
      ContourPlot[
        tKnown[b][x, y],
        {x, xMinCont, xMaxCont}, {y, -yMaxCont, yMaxCont},
        Contours -> numTo1 + numBeyond1,
        ContourShading -> None,
        ContourStyle -> contStyle,
        PlotRange -> {0, 1 + (1 + numBeyond1) / numTo1}
      ],
      (* Straight contour *)
      Graphics @ {straightStyle,
        Line @ {{xStraight, -yMaxCont}, {xStraight, yMaxCont}}
      }
    ]
  , {{b, bInit, bIt}, bMin, bMax}]
]


(* ::Section:: *)
(*Viable domain*)


(* ::Subsection:: *)
(*Interactive visualiser*)


DynamicModule[
 {aInit, aMin, aMax,
  bInit, bMin, bMax,
  xMin, xMax, yMax,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinNon, xMaxNon, yMaxNon,
  xMinCont, xMaxCont, yMaxCont,
  numTo1, numBeyond1
 },
  (* Values of A *)
  aInit = 0.2;
  aMin = 0.01;
  aMax = 2;
  (* Values of B *)
  bInit = 1;
  bMin = 0;
  bMax = 5;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 5/4;
  yMax = 2;
  (* Margin *)
  eps = 0.1;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin - eps;
  xMaxUnphys = xMax + eps;
  yMaxUnphys = yMax + eps;
  (* Plot range for non-viable domain *)
  xMinNon = xMin - eps;
  xMaxNon = xMax + eps;
  yMaxNon = yMax + eps;
  (* Plot range for contours *)
  xMinCont = xMin - eps;
  xMaxCont = xMax + eps;
  yMaxCont = yMax + eps;
  (* Number of contours *)
  numTo1 = 5;
  numBeyond1 = 3;
  (* Plot *)
  Manipulate[
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        ImageSize -> 240,
        PlotLabel -> BoxedLabel[
          Row[
            {aIt == N[a], bIt == N[b]},
            ","
          ]
        ]
      ],
      (* Unphysical domain *)
      RegionPlot[tKnown[b][x, y] < 0,
        {x, xMinUnphys, xMaxUnphys}, {y, -yMaxUnphys, yMaxUnphys},
        BoundaryStyle -> unphysStyle,
        PlotPoints -> 50,
        PlotStyle -> unphysStyle
      ],
      (* Non-viable domain *)
      RegionPlot[vi[a, b][x, y] < 0,
        {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
        BoundaryStyle -> termStyle,
        PlotPoints -> 50,
        PlotStyle -> nonStyle
      ],
      (* Known solution contours *)
      ContourPlot[
        tKnown[b][x, y],
        {x, xMinCont, xMaxCont}, {y, -yMaxCont, yMaxCont},
        Contours -> numTo1 + numBeyond1,
        ContourShading -> None,
        ContourStyle -> contStyle,
        PlotRange -> {0, 1 + (1 + numBeyond1) / numTo1}
      ],
      (* Straight contour *)
      Graphics @ {straightStyle,
        Line @ {{xStraight, -yMaxCont}, {xStraight, yMaxCont}}
      }
    ]
  , {{a, aInit, aIt}, aMin, aMax}
  , {{b, bInit, bIt}, bMin, bMax}]
]


(* Combinations {A, B} which ought to be investigated later *)
{
  {0.2, 1},
  {1, 1},
  {1.2, 1},
  {1, 0.5},
  {0.3, 0.55}
};


(* ::Subsection:: *)
(*Interactive visualiser for \[CapitalPhi] along y = 0*)


DynamicModule[
 {aInit, aMin, aMax,
  bInit, bMin, bMax,
  xMin, xMax
 },
  (* Values of A *)
  aInit = 0.3;
  aMin = 0.01;
  aMax = 3;
  (* Values of B *)
  bInit = 1;
  bMin = 0.1;
  bMax = 3;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 5/4;
  (* Plot *)
  Manipulate[
    Show[
      EmptyFrame[{xMin, xMax}, {-0.1, 1},
        FrameLabel -> {xIt, "\[CapitalPhi]"[yIt == 0]},
        ImageSize -> 360,
        PlotLabel -> BoxedLabel[
          Row[
            {aIt == N[a], bIt == N[b]},
            ","
          ]
        ]
      ],
      Plot[{vi[a, b][x, 0], 0}, {x, xMin, xMax},
        Filling -> {1 -> {2}},
        PlotRange -> Full
      ],
      (* Location of straight contour *)
      Graphics @ {Directive[straightStyle, pointStyle],
        Point @ {xStraight, 0}
      }
    ]
  , {{a, aInit, aIt}, aMin, aMax}
  , {{b, bInit, bIt}, bMin, bMax}]
]


(* ::Subsection:: *)
(*Animation for \[CapitalPhi] along y = 0, simple case (B = 1)*)


Module[
 {aMin, aMax, aStep,
  aValues, b,
  xMin, xMax
 },
  (* Values of A *)
  aMin = 1/10;
  aMax = 2;
  aStep = 1/10;
  aValues = Range[aMin, aMax, aStep];
  (* Value of B *)
  b = 1;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 5/4;
  (* Animation *)
  Table[
    Show[
      EmptyFrame[{xMin, xMax}, {-0.1, 1},
        FrameLabel -> {xIt, "\[CapitalPhi]"[yIt == 0]},
        ImageSize -> 360,
        PlotLabel -> BoxedLabel[
          Row[
            {aIt == N[a], bIt == N[b]},
            ","
          ]
        ]
      ],
      (* \[CapitalPhi] *)
      Plot[{vi[a, b][x, 0], 0}, {x, xMin, xMax},
        Filling -> {1 -> {2}},
        PlotRange -> Full
      ],
      (* Location of straight contour *)
      Graphics @ {Directive[straightStyle, pointStyle],
        Point @ {xStraight, 0}
      },
      (* Non-trivial critical terminal point x_0 *)
      Graphics @ {Directive[Red, pointStyle = PointSize[Large]],
        Point @ {x0Simp[a], 0}
      }
    ]
  , {a, aValues}]
] // Ex["cosine_simple-phi.gif", gifOpts]


(* ::Subsection:: *)
(*Animation for viable domain, simple case (B = 1)*)


Module[
 {aMin, aMax, aStep,
  aValues, b,
  xMin, xMax, yMax,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinNon, xMaxNon, yMaxNon,
  xMinCont, xMaxCont, yMaxCont,
  unphysicalDomain,
  numTo1, numBeyond1,
  genericContours, straightContour
 },
  (* Values of A *)
  aMin = 1/10;
  aMax = 2;
  aStep = 1/10;
  aValues = Range[aMin, aMax, aStep];
  (* Value of B *)
  b = 1;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 5/4;
  yMax = 2;
  (* Margin *)
  eps = 0.15;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin - eps;
  xMaxUnphys = xMax + eps;
  yMaxUnphys = yMax + eps;
  (* Plot range for non-viable domain *)
  xMinNon = xMin - eps;
  xMaxNon = xMax + eps;
  yMaxNon = yMax + eps;
  (* Plot range for contours *)
  xMinCont = xMin - eps;
  xMaxCont = xMax + eps;
  yMaxCont = yMax + eps;
  (*
    NOTE: since B is fixed at 1,
    the unphysical domain and known solution contours
    do not need to be redrawn for different values of A.
   *)
  (* Unphysical domain *)
  unphysicalDomain =
    RegionPlot[tKnown[b][x, y] < 0,
      {x, xMinUnphys, xMaxUnphys}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> unphysStyle,
      PlotPoints -> 50,
      PlotStyle -> unphysStyle
    ];
  (* Known solution contours *)
  numTo1 = 5;
  numBeyond1 = 3;
  genericContours =
    ContourPlot[
      tKnown[b][x, y],
      {x, xMinCont, xMaxCont}, {y, -yMaxCont, yMaxCont},
      Contours -> numTo1 + numBeyond1,
      ContourShading -> None,
      ContourStyle -> contStyle,
      PlotRange -> {0, 1 + (1 + numBeyond1) / numTo1}
    ];
  straightContour =
    Graphics @ {straightStyle,
      Line @ {{xStraight, -yMaxCont}, {xStraight, yMaxCont}}
    };
  (* Animation *)
  Table[
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        ImageSize -> 180,
        PlotLabel -> BoxedLabel[
          Row[
            {aIt == N[a], bIt == N[b]},
            ","
          ]
        ]
      ],
      (* Unphysical domain *)
      unphysicalDomain,
      (* Non-viable domain *)
      RegionPlot[vi[a, b][x, y] < 0,
        {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
        BoundaryStyle -> termStyle,
        PlotPoints -> 50,
        PlotStyle -> nonStyle
      ],
      (* Known solution contours (generic) *)
      genericContours,
      straightContour,
      (* Critical terminal point (x_0, 0) *)
      Graphics @ {Directive[simp0Style, pointStyle],
        Point @ {x0Simp[a], 0}
      },
      (* Known solution contour through (x_0, 0) *)
      ContourPlot[
        tKnown[b][x, y] == tKnown[b][x0Simp[a], 0],
        {x, xMinCont, xMaxCont}, {y, -yMaxCont, yMaxCont},
        ContourStyle -> simp0Style,
        PlotRange -> {0, 1}
      ]
    ]
  , {a, aValues}]
] // Ex["cosine_simple-viable.gif", gifOpts]


(* ::Section:: *)
(*Simple case (B = 1) plots*)


(* ::Subsection:: *)
(*Starting points for boundary tracing*)


Module[
 {b,
  xMin, xMax, yMax,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinNon, xMaxNon, yMaxNon,
  xMinCont, xMaxCont, yMaxCont,
  unphysicalDomain,
  numTo1, numBeyond1, genericContours,
  xStraight, straightContour,
  idList
 },
  (* Value of B *)
  b = 1;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 5/4;
  yMax = 2;
  (* Margin *)
  eps = 0.1;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin - eps;
  xMaxUnphys = xMax + eps;
  yMaxUnphys = yMax + eps;
  (* Plot range for non-viable domain *)
  xMinNon = xMin - eps;
  xMaxNon = xMax + eps;
  yMaxNon = yMax + eps;
  (* Plot range for contours *)
  xMinCont = xMin - eps;
  xMaxCont = xMax + eps;
  yMaxCont = yMax + eps;
  (*
    NOTE: since B is fixed at 1,
    the unphysical domain and known solution contours
    do not need to be redrawn for different values of A.
   *)
  (* Unphysical domain *)
  unphysicalDomain =
    RegionPlot[tKnown[b][x, y] < 0,
      {x, xMinUnphys, xMaxUnphys}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> unphysStyle,
      PlotPoints -> 50,
      PlotStyle -> unphysStyle
    ];
  (* Known solution contours (generic & straight) *)
  numTo1 = 5;
  numBeyond1 = 3;
  genericContours =
    ContourPlot[
      tKnown[b][x, y],
      {x, xMinCont, xMaxCont}, {y, -yMaxCont, yMaxCont},
      Contours -> numTo1 + numBeyond1,
      ContourShading -> None,
      ContourStyle -> contStyle,
      PlotRange -> {0, 1 + (1 + numBeyond1) / numTo1}
    ];
  xStraight = Pi/2;
  straightContour = Graphics @ {straightStyle,
    Line @ {{xStraight, -yMaxCont}, {xStraight, yMaxCont}}
  };
  (* Group names *)
  idList = {"contour", "hyperbolic"};
  (* Plots for various A *)
  Table[
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        ImageSize -> 240,
        PlotLabel -> BoxedLabel[
          Row[
            {aIt == N[a], bIt == N[b]},
            ","
          ]
        ]
      ],
      (* Unphysical domain *)
      unphysicalDomain,
      (* Non-viable domain *)
      RegionPlot[vi[a, b][x, y] < 0,
        {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
        BoundaryStyle -> termStyle,
        PlotPoints -> 50,
        PlotStyle -> nonStyle
      ],
      (* Known solution contours (generic) *)
      genericContours,
      straightContour,
      (* Starting points *)
      ListPlot[
        Table[startXYSimp[a][id], {id, idList}],
        LabelingFunction -> Function @ Placed[
          #2[[2]],
          Center
        ],
        PlotLegends -> idList,
        PlotStyle -> Directive[Opacity[0.7], pointStyle]
      ]
    ]
    // Ex @ FString[
      "cosine_simple-a-{aNamesSimp[a]}-traced-starting.pdf"
    ]
  , {a, aValuesSimp}]
]


(* ::Subsection:: *)
(*Traced boundaries*)


Module[
 {b,
  xMin, xMax, yMax,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinNon, xMaxNon, yMaxNon,
  unphysicalDomain,
  xStraight, straightContour,
  idList
 },
  (* Value of B *)
  b = 1;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 5/4;
  yMax = 2;
  (* Margin *)
  eps = 0.1;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin - eps;
  xMaxUnphys = xMax + eps;
  yMaxUnphys = yMax + eps;
  (* Plot range for non-viable domain *)
  xMinNon = xMin - eps;
  xMaxNon = xMax + eps;
  yMaxNon = yMax + eps;
  (*
    NOTE: since B is fixed at 1,
    the unphysical domain and known solution contours
    do not need to be redrawn for different values of A.
   *)
  (* Unphysical domain *)
  unphysicalDomain =
    RegionPlot[tKnown[b][x, y] < 0,
      {x, xMinUnphys, xMaxUnphys}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> unphysStyle,
      PlotPoints -> 50,
      PlotStyle -> unphysStyle
    ];
  (* Straight contour *)
  xStraight = Pi/2;
  straightContour = Graphics @ {straightStyle,
    Line @ {{xStraight, -yMaxNon}, {xStraight, yMaxNon}}
  };
  (* Group names *)
  idList = {"contour", "hyperbolic"};
  (* Plots for various A *)
  Table[
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        ImageSize -> 240,
        PlotLabel -> BoxedLabel[
          Row[
            {aIt == N[a], bIt == N[b]},
            ","
          ]
        ]
      ],
      (* Unphysical domain *)
      unphysicalDomain,
      (* Non-viable domain *)
      RegionPlot[vi[a, b][x, y] < 0,
        {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
        BoundaryStyle -> termStyle,
        PlotPoints -> 50,
        PlotStyle -> nonStyle
      ],
      (* Straight contour *)
      straightContour,
      (* Traced boundaries (general) *)
      Table[
        Table[
          ParametricPlot[
            xy[s]
              // Through
              // {#, {#[[1]], -#[[2]]}} &
              // Evaluate,
            {s, DomainStart[xy], DomainEnd[xy]},
            PlotStyle -> {upperStyle, lowerStyle}
          ]
        , {xy, xyTraSimp[a][id]}]
      , {id, idList}],
      (* Traced boundaries (hyperbolic) *)
      Table[
        Table[
          ParametricPlot[
            xy[s]
              // Through
              // {#, {#[[1]], -#[[2]]}} &
              // Evaluate,
            {s, DomainStart[xy], DomainEnd[xy]},
            PlotStyle -> glowStyle
          ]
        , {xy, xyTraSimp[a][id]}]
      , {id, {"hyperbolic"}}]
    ]
    // Ex @ FString[
      "cosine_simple-a-{aNamesSimp[a]}-traced-full.pdf"
    ]
  , {a, aValuesSimp}]
]


(* ::Text:: *)
(*The terminal curve is practically a traced boundary.*)
(*For each A < 1 there will be some (x_i, y_i) where it inflects;*)
(*if x_i > \[Pi]/2, then we will have a convex domain*)
(*bounded by that traced boundary and by x = pi/2,*)
(*with thermal radiation and constant temperature along them respectively.*)
(*For a (contrived) practical situation, imagine rod at constant temperature*)
(*with a knob affixed at the end, and radiating into vacuum.*)


(* ::Subsection:: *)
(*Check boundary condition along terminal curve*)


With[{x = \[FormalX], y = \[FormalY], a = \[FormalCapitalA]},
  Module[
   {aValues, b, t, gradT,
    n, fluxL, fluxR, res,
    sMax,
    xTerm, yTerm
   },
    (* Values of A *)
    aValues = {1/10^3, 1/100, 1/10, 1/2, 3/4, 1, 5/4, 3/2};
    (* Value of B *)
    b = 1;
    (* Known solution T *)
    t = tKnown[b][x, y];
    (* Gradient of T *)
    gradT = Grad[t, {x, y}];
    (* Normal vector n *)
    n = -{y', -x'};
    (* Left hand side of boundary condition *)
    fluxL = n . gradT;
    (* Right hand side of boundary condition *)
    fluxR = -t^4 / a;
    (* Residual of boundary condition *)
    res["abs"] = fluxL - fluxR;
    res["rel"] = fluxR / fluxL - 1;
    (* Plots of residuals *)
    sMax = 5;
    Table[
      Plot[
        Table[
          (* Terminal curve x == x(s), y == y(s) *)
          {xTerm, yTerm} = xyTermSimp[a];
          (* Evaluate residual therealong *)
          res[type] /. {
            x' -> xTerm'[s],
            y' -> yTerm'[s],
            x -> xTerm[s],
            y -> yTerm[s]
          }
        , {a, aValues}] // Evaluate,
        {s, -sMax, sMax},
        AxesLabel -> {sIt, FString @ "Residual ({type})"},
        ImageSize -> 360,
        PlotLabel -> "Along terminal curve:",
        PlotLegends -> LineLegend[aValues // N, LegendLabel -> aIt],
        PlotOptions[Axes] // Evaluate
      ] // Ex @ FString[
        "cosine_simple-terminal-boundary-residual-{type}.pdf"
      ]
    , {type, {"abs", "rel"}}]
  ]
]


(* ::Subsection:: *)
(*Check boundary condition along candidate traced boundary*)


(* ::Text:: *)
(*(The candidate traced boundary is the traced boundary through (x_0, 0),*)
(*which is almost the same as the terminal curve.)*)


With[{x = \[FormalX], y = \[FormalY], a = \[FormalCapitalA]},
  Module[
   {aValues, b, t, gradT,
    n, fluxL, fluxR, res,
    sMax,
    xTerm, yTerm
   },
    (* Values of A *)
    aValues = {1/10^3, 1/100, 1/10, 1/2, 3/4, 1, 5/4, 3/2};
    (* Value of B *)
    b = 1;
    (* Known solution T *)
    t = tKnown[b][x, y];
    (* Gradient of T *)
    gradT = Grad[t, {x, y}];
    (* Normal vector n *)
    n = {y', -x'};
    (* Left hand side of boundary condition *)
    fluxL = n . gradT;
    (* Right hand side of boundary condition *)
    fluxR = -t^4 / a;
    (* Residual of boundary condition *)
    res["abs"] = fluxL - fluxR;
    res["rel"] = fluxR / fluxL - 1;
    (* Plots of residuals *)
    sMax = 5;
    Table[
      Plot[
        Table[
          (* Candidate boundary x == x(s), y == y(s) *)
          {xTerm, yTerm} = xyTraCandSimp[a];
          (* Evaluate residual therealong *)
          res[type] /. {
            x' -> xTerm'[s],
            y' -> yTerm'[s],
            x -> xTerm[s],
            y -> yTerm[s]
          } /. {
            s -> Abs[s]
          }
        , {a, aValues}] // Evaluate,
        {s, -sMax, sMax},
        AxesLabel -> {sIt, FString @ "Residual ({type})"},
        ImageSize -> 360,
        PlotLabel -> "Along terminal curve:",
        PlotLegends -> LineLegend[aValues // N, LegendLabel -> aIt],
        PlotRange -> Full,
        PlotOptions[Axes] // Evaluate
      ] // Ex @ FString[
        "cosine_simple-candidate-boundary-residual-{type}.pdf"
      ]
    , {type, {"abs", "rel"}}]
  ]
]


(* ::Subsection:: *)
(*Curvature algebra*)


(* ::Text:: *)
(*See (r5.21) (Page r5-5).*)


With[{y = \[FormalY], a = \[FormalCapitalA]},
  Block[{$Assumptions = 0 < a <= 1 && 0 < y < ArcCosh[1/a]},
    curTra[a, 1][Pi/2, y] == (
      a^2 Cosh[y] / Sqrt[a^2 Cosh[y]^2 - 1] * (
        2 Sinh[y] - Cosh[y]^2 (a^2 Sinh[y] + 4 Sqrt[a^2 Cosh[y]^2 - 1])
      )
    )
    // FullSimplify
  ]
]


(* ::Text:: *)
(*See (r5.22) & (r5.23) (Page r5-6).*)


With[{y = \[FormalY], a = \[FormalCapitalA], s = \[FormalCapitalS]},
  2 Sinh[y] - Cosh[y]^2 (a^2 Sinh[y] + 4 Sqrt[a^2 Cosh[y]^2 - 1]) == (
    2 s - (1 + s^2) (a^2 s + 4 Sqrt[a^2 (1 + s^2) - 1])
  )
    /. {s -> Sinh[y]}
    // FullSimplify
]


(* ::Subsection:: *)
(*Horizontal coordinate x at critical y = y(A)*)


(* ::Text:: *)
(*A plot of x(y(A)) - \[Pi]/2, whose zero A = A_i was sought.*)
(*See (r5.25) (Page r5-6).*)


Module[
 {aMin, aMax, aInfl
 },
  (* Plot range *)
  aMin = 0.15;
  aMax = 1;
  (* Value of A_i *)
  aInfl = aInflSimp;
  (* Plot *)
  Show[
    (* A vs x(y(A)) - \[Pi]/2 *)
    Plot[
      xTraCandSimp[a] @ yCurCritSimp[a] - Pi / 2 // Evaluate,
      {a, aMin, aMax},
      AxesLabel -> {aIt, xIt @ yIt[aIt] - Pi / 2},
      ImageSize -> 360,
      PlotOptions[Axes] // Evaluate
    ],
    (* A_i *)
    Graphics @ {inflDotStyle,
      Point @ {aInfl, 0}
    },
    Graphics @ Text[
      Subscript[aIt, "i"] == N[aInfl]
        // textStyle,
      {aInfl, 0},
      {0, -1.3}
    ]
  ]
] // Ex["cosine_simple-candidate-x-minus-half-pi.pdf"]


(* ::Subsection:: *)
(*Approximate A_i (more curvature algebra)*)


(* ::Text:: *)
(*Check (r5.28) (Page r5-6).*)


With[{a = \[FormalCapitalA]},
  vi[a, 1][Pi/2, ArcSech[a]] == 0
]


(* ::Subsubsection:: *)
(*Curvature of \[CapitalPhi] = 0 at (x, y) = (\[Pi]/2, arcsech(A))*)


(* ::Text:: *)
(*Recall that the terminal curve \[CapitalPhi] = 0*)
(*is almost indistinguishable from the candidate boundary.*)
(*See (r5.29) (Page r5-7).*)


With[{x = \[FormalX], y = \[FormalY], a = \[FormalCapitalA]},
  Module[{div, grad, cur},
    (* Abbreviations *)
    div = Div[#, {x, y}] &;
    grad = Grad[#, {x, y}] &;
    (* Curvature *)
    cur = div @ Normalize @ grad @ vi[a, 1][x, y];
    (* Evaluate at crossing *)
    cur /. {x -> Pi/2, y -> ArcSech[a]}
      // FullSimplify[#, 0 < a < 1] &
  ]
]


Module[
 {aMin, aMax, aInfl
 },
  (* Plot range *)
  aMin = 0;
  aMax = 1;
  (* Value of A_i(approx) *)
  aInfl = aInflSimpApprox;
  (* Plot *)
  Show[
    (* A vs \[Kappa] *)
    Plot[
      a (a^6 - a^4 + 44 a^2 - 28) / (16 + a^2 - a^4)^(3/2),
      {a, aMin, aMax},
      AxesLabel -> {aIt, "\[Kappa]"},
      ImageSize -> 360,
      PlotLabel -> Column @ {
        "Curvature of terminal curve",
        Row @ {"at ", ""[xIt == Pi/2, yIt == ArcSech[aIt]]}
      },
      PlotOptions[Axes] // Evaluate
    ],
    (* A_i *)
    Graphics @ {inflDotStyle,
      Point @ {aInfl, 0}
    },
    Graphics @ Text[
      Subscript[aIt, "i"] == N[aInfl]
        // textStyle,
      {aInfl, 0},
      {0, -1.3}
    ]
  ]
] // Ex["cosine_simple-terminal-curvature-x-at-half-pi.pdf"]


(* ::Subsection:: *)
(*Convex domains without known solution*)


Module[
 {b,
  xRad, yEnd,
  xMin, xMax, yMax
 },
  Table[
    (* Value of B *)
    b = 1;
    (* Radiation boundary x == x(y) for convex domain *)
    xRad = xTraCandSimp[a, True];
    yEnd = DomainEnd[xRad];
    xRad = Function[{y}, xRad[Abs @ y] // Evaluate];
    (* Plot range *)
    xMin = xStraight - 0.6 yEnd;
    xMax = xStraight + 0.6 yEnd;
    yMax = yEnd;
    (* Plot *)
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        ImageSize -> 240,
        PlotLabel -> BoxedLabel[aIt == N[a]]
      ],
      (* Convex domain *)
      ParametricPlot[
        {
          {xRad[y], y},
          {xStraight, y}
        }, {y, -yEnd, yEnd},
        PlotStyle -> convexStyle
      ]
    ] // Ex @ FString[
      "cosine_simple-traced-{aNamesSimpConvex[a]}.pdf"
    ]
  , {a, aValuesSimpConvex}]
]


(* ::Subsection:: *)
(*Convex domain aspect ratio*)


Module[{aMin, aMax},
  (* Plot range *)
  aMin = aInflSimp;
  aMax = 1;
  (* Plot *)
  Show[
    (* A vs \[Kappa] *)
    Plot[
      aspectRatioSimp[a],
      {a, aMin, aMax},
      AxesLabel -> {aIt, Null},
      ImageSize -> 360,
      PlotLabel -> "Aspect ratio of convex domains",
      PlotRange -> {0, Automatic},
      PlotOptions[Axes] // Evaluate
    ]
  ]
] // Ex["cosine_simple-traced-convex-aspect-ratio.pdf"]
