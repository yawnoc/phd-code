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
(*Italicised symbols*)


aIt = Italicised["A"];
bIt = Italicised["B"];


(* ::Subsection:: *)
(*Global styles for plots*)


nonStyle = Directive[Opacity[0.7], LightGray];
termStyle = Pink;
unphysStyle = Black;


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
  xMin, xMax, yMax
 },
  (* Values of B *)
  bInit = 1;
  bMin = 0;
  bMax = 5;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2;
  yMax = 2;
  (* Plot *)
  Manipulate[
    RegionPlot[
      tKnown[b][x, y] < 0,
      {x, xMin, xMax}, {y, -yMax, yMax},
      BoundaryStyle -> unphysStyle,
      PlotLabel -> BoxedLabel[bIt == N[b]],
      PlotOptions[Frame] // Evaluate,
      PlotStyle -> unphysStyle
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
  xMin, xMax, yMax
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
  (* Plot *)
  Manipulate[
    Show[
      EmptyFrame[{xMin, xMax}, {-xMax, xMax},
        PlotLabel -> BoxedLabel[
          Row[
            {aIt == N[a], bIt == N[b]},
            ","
          ]
        ]
      ],
      (* Unphysical domain *)
      RegionPlot[tKnown[b][x, y] < 0,
        {x, xMin, xMax}, {y, -yMax, yMax},
        BoundaryStyle -> unphysStyle,
        PlotStyle -> unphysStyle
      ],
      (* Non-viable domain *)
      RegionPlot[vi[a, b][x, y] < 0,
        {x, xMin, xMax}, {y, -yMax, yMax},
        BoundaryStyle -> termStyle,
        PlotStyle -> nonStyle
      ]
    ]
  , {{a, aInit, aIt}, aMin, aMax}
  , {{b, bInit, bIt}, bMin, bMax}]
]
