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


tKnown[b_][x_, y_] := 1 - b Cos[x] Cosh[y];


(* ::Subsubsection:: *)
(*P = \[PartialD]T/\[PartialD]x, Q = \[PartialD]T/\[PartialD]y*)


p[b_][x_, y_] := D[tKnown[b][x, y], x] // Evaluate;
q[b_][x_, y_] := D[tKnown[b][x, y], y] // Evaluate;


(* ::Subsection:: *)
(*Flux function F*)


f[a_, b_][x_, y_] := -tKnown[b][x, y]^4 / a // Evaluate;


(* ::Subsection:: *)
(*Viability \[CapitalPhi]*)


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


(* ::Subsection:: *)
(*Italicised symbols*)


aIt = Italicised["A"];
bIt = Italicised["B"];
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


contStyle = LightGray;
nonStyle = Directive[Opacity[0.7], LightGray];
termStyle = Pink;
unphysStyle = Black;


simp0Style = Blue;


pointStyle = PointSize[Large];


(* ::Section:: *)
(*Known solution*)


(* ::Subsection:: *)
(*Algebra check for \[CapitalPhi]*)


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
  xMinCont, xMaxCont, yMaxCont
 },
  (* Values of B *)
  bInit = 1;
  bMin = 0;
  bMax = 5;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2;
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
        ContourShading -> None,
        ContourStyle -> contStyle,
        PlotRange -> {0, 1}
      ]
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
  xMinCont, xMaxCont, yMaxCont
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
  xMax = Pi/2;
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
        ContourShading -> None,
        ContourStyle -> contStyle,
        PlotRange -> {0, 1}
      ]
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
  xMax = Pi/2;
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
      ]
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
  xMax = Pi/2;
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
  genericContours
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
  xMax = Pi/2;
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
  (* Known solution contours (generic) *)
  genericContours =
    ContourPlot[
      tKnown[b][x, y],
      {x, xMinCont, xMaxCont}, {y, -yMaxCont, yMaxCont},
      ContourShading -> None,
      ContourStyle -> contStyle,
      PlotRange -> {0, 1}
    ];
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
