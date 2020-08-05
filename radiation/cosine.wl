(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< NDSolve`FEM`
<< Conway`
<< FigureStyles`
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
      /. Last @ Solve[
        {
          (1 - c^2) - (1 - c)^8 / a^2 == 0,
          c < 1
        }, c, Reals
      ]
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
(*General case (B arbitrary)*)


(* ::Subsubsection:: *)
(*Critical terminal points along y = 0*)


(* ::Text:: *)
(*Observe that \[CapitalPhi](y = 0) = B^2 (1 - c^2) - (1 - B c)^8 / A^2, where c = cos(x).*)
(*For a given A, there exists a constant B_\[Natural](A) such that*)
(*only for B > B_\[Natural](A) will there be two critical terminal x-values,*)
(*x_\[Flat] < x_\[Sharp], otherwise there will be none.*)


(* ::Subsubsubsection:: *)
(*Exploratory analysis: solutions of the polynomial*)


With[{c = \[FormalC], a = \[FormalCapitalA], b = \[FormalCapitalB]},
  Solve[
    {
      b^2 (1 - c^2) - (1 - b c)^8 / a^2 == 0,
      a > 0,
      b > 0
    }, c, Reals
  ]
]


(* ::Subsubsubsection:: *)
(*Critical value B_\[Natural](A)*)


(* ::Text:: *)
(*From the above output,*)
(*we see that the critical constant B_\[Natural](A) is given by*)
(*the edge of the inequality A > Root[-65536 + ... + 729 B^2 #1^4 &, 2].*)
(*Inverting that expression for B in terms of A:*)


With[{a = \[FormalCapitalA], b = \[FormalCapitalB]},
  Solve[
    {
      a == Root[
        -65536 + 458752 b^2 - 1376256 b^4 + 2293760 b^6 - 2293760 b^8
        + 1376256 b^10 - 458752 b^12 + 65536 b^14
        + (-729 + 77200 b^2 + 449856 b^4 + 283392 b^6 + 13824 b^8) #1^2
        + 729 b^2 #1^4 &,
        2
      ],
      b > 0
    }, b, Reals
  ]
]


(* ::Text:: *)
(*Therefore:*)


bNat[a_] := Root[
  -65536 - 729 a^2
  + (458752 + 77200 a^2 + 729 a^4) #1^2
  + (-1376256 + 449856 a^2) #1^4
  + (2293760 + 283392 a^2) #1^6
  + (-2293760 + 13824 a^2) #1^8
  + 1376256 #1^10
  - 458752 #1^12
  + 65536 #1^14 &,
  2
]


(* ::Subsubsubsection:: *)
(*Critical terminal points x_\[Flat] and x_\[Sharp]*)


(* ::Text:: *)
(*From the above output (see Exploratory analysis: solutions of the polynomial),*)
(*we see that the two values of critical terminal cos(x)*)
(*are given by the Root[1 - A^2 B^2 - ... + B^8 #1^8 &, ...] expressions:*)


polyC[a_, b_] := Function[{c},
  1 - a^2 b^2
  - 8 b c
  + (28 b^2 + a^2 b^2) c^2
  - 56 b^3 c^3
  + 70 b^4 c^4
  - 56 b^5 c^5
  + 28 b^6 c^6
  - 8 b^7 c^7
  + b^8 c^8
];


xFlat[a_, b_] := ArcCos @ Root[polyC[a, b], 2] // Evaluate;


xSharp[a_, b_] := ArcCos @ Root[polyC[a, b], 1] // Evaluate;


(* ::Subsubsubsection:: *)
(*Critical x_\[Natural] at B = B_\[Natural](A)*)


xNat[a_] := xSharp[a, bNat[a]] // Evaluate;


(* ::Subsubsection:: *)
(*Representative values of A*)


aValuesGen = {3, 10};


aNamesGen = AssociationThread[
  aValuesGen,
  {"three", "ten"}
];


(* ::Subsubsection:: *)
(*Representative values of B for a given A*)


(* ::Text:: *)
(*The names of the regimes are:*)
(*(1) gentle, B < B_\[Natural](A)*)
(*(2) gentle-to-fair, B = B_\[Natural](A)*)
(*(3) fair, B_\[Natural](A) < B < 1*)
(*(4) fair-to-steep (or simple), B = 1*)
(*(5) steep, B > 1*)
(*See Page r5-8.*)


(* Gentle regime: choose B == 9/10 B_\[Natural](A) *)
bValueGen["gentle"][a_] := 9/10 bNat[a];


(* Gentle-to-fair transition: B == B_\[Natural](A) *)
bValueGen["gentle_fair"][a_] := bNat[a];


(* Fair regime: B_\[Natural] < B < 1 such that x_\[Flat] == 3/4 x_\[Natural] *)
bValueGen["fair"][a_] :=
  SeekRootBisection[
    xFlat[a, #] - 3/4 xNat[a] &,
    {bNat[a], 1}
  ];


(* Fair-to-steep transition: B == 1 *)
bValueGen["fair_steep"][a_] := 1;


(* Steep regime: choose B == 2 *)
bValueGen["steep"][a_] := 2;


(* ::Subsubsection:: *)
(*Representative pair (A, B) for an asymmetric domain*)


aAsymm = 12;
bAsymm = 105/100 bNat[aAsymm];


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
        (* Contour *)
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
        (* Starting points along the contour *)
        Table[
          xyContour[s] // Through // Rationalize[#, 0] &
        , {s, Subdivide[-sMax, sMax, nMax]}]
      ]
    ];
, {a, aValuesSimp}];


(* ::Subsubsection:: *)
(*General case (B arbitrary)*)


(* ::Subsubsubsection:: *)
(*Gentle regime B < B_\[Natural](A)*)


Table[
  Module[
   {yReflect,
    regime,
    b, yMax,
    tValue, xInit, yInit,
    sMax, xyContour,
    sStart, sEnd,
    nMax
   },
    (* Reflect y coordinate *)
    yReflect = # * {1, -1} &;
    (* Regime name *)
    regime = "gentle";
    (* Value of B *)
    b = bValueGen[regime][a];
    (* Range for y *)
    yMax = 4;
    (*
      ----------------------------------------------------------------
      Starting points along (probably) disconnected contour
        T = 2/3 T(x = x_\[Natural], y = 0)
      (For the current values A = 3, 10, this contour is disconnected,
      but I say "probably" since there might be some extreme choice of A & B
      for which the contour is not disconnected.)
      ----------------------------------------------------------------
    *)
    (* Value of T along contour *)
    tValue = 2/3 tKnown[b][xNat[a], 0];
    (* Seek starting point (x == 1, y) therealong *)
    xInit = 1;
    yInit = SeekRoot[
      tKnown[b][xInit, #] - tValue &,
      {0, yMax}
    ];
    (* (Probable) upper bound for arc length traversed *)
    sMax = 6;
    (* Contour (y > 0 half) *)
    xyContour =
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyContourSystem[b],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {
                tKnown[b][x[s], y[s]] < 0,
                x[s] < 0,
                y[s] > yMax
              },
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax},
          NoExtrapolation
        ]
      ];
    (* Actual arc length traversed *)
    sStart = DomainStart[xyContour];
    sEnd = DomainEnd[xyContour];
    (* Return starting points *)
    nMax = 6;
    startXYGen[a][regime]["disconnected"] = (
      Table[
        xyContour[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[sStart, sEnd, nMax]}]
        // Rest
        // Most
        // Join[#, yReflect /@ #] &
    );
    (*
      ----------------------------------------------------------------
      Starting points along (probably) connected contour
        T = T(x = x_\[Natural], y = 0)
      ----------------------------------------------------------------
    *)
    (* Starting point (x == x_\[Natural], y == 0) *)
    xInit = xNat[a];
    yInit = 0;
    (* (Probable) upper bound for arc length traversed *)
    sMax = 6;
    (* Contour *)
    xyContour =
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyContourSystem[b],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {Abs @ y[s] > yMax},
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax},
          NoExtrapolation
        ]
      ];
    (* Actual arc length traversed *)
    sStart = DomainStart[xyContour];
    sEnd = DomainEnd[xyContour];
    (* Return starting points within the viable domain *)
    nMax = 12;
    startXYGen[a][regime]["connected"] = (
      Table[
        xyContour[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[sStart, sEnd, nMax]}]
        // Rest
        // Most
        // Cases[{x_, y_} /; vi[a, b][x, y] > 0]
    );
  ]
, {a, aValuesGen}];


(* ::Subsubsubsection:: *)
(*Gentle-to-fair transition B = B_\[Natural](A)*)


Table[
  Module[
   {yReflect,
    regime,
    b, yMax,
    tValue, xInit, yInit,
    sMax, xyContour,
    sStart, sEnd,
    nMax
   },
    (* Reflect y coordinate *)
    yReflect = # * {1, -1} &;
    (* Regime name *)
    regime = "gentle_fair";
    (* Value of B *)
    b = bValueGen[regime][a];
    (* Range for y *)
    yMax = 4;
    (*
      ----------------------------------------------------------------
      Starting points along (probably) disconnected contour
        T = 2/3 T(x = x_\[Natural], y = 0)
      (For the current values A = 3, 10, this contour is disconnected,
      but I say "probably" since there might be some extreme choice of A & B
      for which the contour is not disconnected.)
      ----------------------------------------------------------------
    *)
    (* Value of T along contour *)
    tValue = 2/3 tKnown[b][xNat[a], 0];
    (* Seek starting point (x == 1, y) therealong *)
    xInit = 1;
    yInit = SeekRoot[
      tKnown[b][xInit, #] - tValue &,
      {0, yMax}
    ];
    (* (Probable) upper bound for arc length traversed *)
    sMax = 6;
    (* Contour (y > 0 half) *)
    xyContour =
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyContourSystem[b],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {
                tKnown[b][x[s], y[s]] < 0,
                x[s] < 0,
                y[s] > yMax
              },
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax},
          NoExtrapolation
        ]
      ];
    (* Actual arc length traversed *)
    sStart = DomainStart[xyContour];
    sEnd = DomainEnd[xyContour];
    (* Return starting points *)
    nMax = 6;
    startXYGen[a][regime]["disconnected"] = (
      Table[
        xyContour[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[sStart, sEnd, nMax]}]
        // Rest
        // Most
        // Join[#, yReflect /@ #] &
    );
    (*
      ----------------------------------------------------------------
      Starting points along connected contour
        T = T(x = x_\[Natural], y = 0)
      ----------------------------------------------------------------
    *)
    (* Starting point (x == x_\[Natural], y == 0) *)
    xInit = xNat[a];
    yInit = 0;
    (* (Probable) upper bound for arc length traversed *)
    sMax = 6;
    (* Contour *)
    xyContour =
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyContourSystem[b],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {Abs @ y[s] > yMax},
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax},
          NoExtrapolation
        ]
      ];
    (* Actual arc length traversed *)
    sStart = DomainStart[xyContour];
    sEnd = DomainEnd[xyContour];
    (* Return starting points within the viable domain *)
    nMax = 12;
    startXYGen[a][regime]["connected"] = (
      Table[
        xyContour[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[sStart, sEnd, nMax]}]
        // Rest
        // Most
    );
    (*
      ----------------------------------------------------------------
      Hyperbolic critical terminal point (x = x_\[Natural], y = 0)
      ----------------------------------------------------------------
    *)
    startXYGen[a][regime]["hyperbolic"] = {
      {xNat[a], 0}
    }
  ]
, {a, aValuesGen}];


(* ::Subsubsubsection:: *)
(*Fair regime B_\[Natural](A) < B < 1*)


Table[
  Module[
   {yReflect,
    regime,
    b, yMax,
    tValue, xInit, yInit,
    sMax, xyContour,
    sStart, sEnd,
    nMax
   },
    (* Reflect y coordinate *)
    yReflect = # * {1, -1} &;
    (* Regime name *)
    regime = "fair";
    (* Value of B *)
    b = bValueGen[regime][a];
    (* Range for y *)
    yMax = 4;
    (*
      ----------------------------------------------------------------
      Starting points along (probably) disconnected contour
        T = 2/3 T(x = x_\[Natural], y = 0)
      (For the current values A = 3, 10, this contour is disconnected,
      but I say "probably" since there might be some extreme choice of A & B
      for which the contour is not disconnected.)
      ----------------------------------------------------------------
    *)
    (* Value of T along contour *)
    tValue = 2/3 tKnown[b][xNat[a], 0];
    (* Seek starting point (x == 1, y) therealong *)
    xInit = 1;
    yInit = SeekRoot[
      tKnown[b][xInit, #] - tValue &,
      {0, yMax}
    ];
    (* (Probable) upper bound for arc length traversed *)
    sMax = 6;
    (* Contour (y > 0 half) *)
    xyContour =
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyContourSystem[b],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {
                tKnown[b][x[s], y[s]] < 0,
                x[s] < 0,
                y[s] > yMax
              },
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax},
          NoExtrapolation
        ]
      ];
    (* Actual arc length traversed *)
    sStart = DomainStart[xyContour];
    sEnd = DomainEnd[xyContour];
    (* Return starting points *)
    nMax = 6;
    startXYGen[a][regime]["disconnected"] = (
      Table[
        xyContour[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[sStart, sEnd, nMax]}]
        // Rest
        // Most
        // Join[#, yReflect /@ #] &
    );
    (*
      ----------------------------------------------------------------
      Starting points along (probably) connected contour
        T = T(x = x_\[Natural], y = 0)
      ----------------------------------------------------------------
    *)
    (* Starting point (x == x_\[Natural], y == 0) *)
    xInit = xNat[a];
    yInit = 0;
    (* (Probable) upper bound for arc length traversed *)
    sMax = 6;
    (* Contour *)
    xyContour =
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyContourSystem[b],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {Abs @ y[s] > yMax},
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax},
          NoExtrapolation
        ]
      ];
    (* Actual arc length traversed *)
    sStart = DomainStart[xyContour];
    sEnd = DomainEnd[xyContour];
    (* Return starting points within the viable domain *)
    nMax = 12;
    startXYGen[a][regime]["connected"] = (
      Table[
        xyContour[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[sStart, sEnd, nMax]}]
        // Rest
        // Most
    );
    (*
      ----------------------------------------------------------------
      Hyperbolic critical terminal points
        (x = x_\[Flat], y = 0)
        (x = x_\[Sharp], y = 0)
      ----------------------------------------------------------------
    *)
    startXYGen[a][regime]["hyperbolic"] = {
      {xFlat[a, b], 0},
      {xSharp[a, b], 0}
    }
  ]
, {a, aValuesGen}];


(* ::Subsubsubsection:: *)
(*Fair-to-steep transition B = 1*)


Table[
  Module[
   {yReflect,
    regime,
    b, yMax,
    tValue, xInit, yInit,
    sMax, xyContour,
    sStart, sEnd,
    nMax
   },
    (* Reflect y coordinate *)
    yReflect = # * {1, -1} &;
    (* Regime name *)
    regime = "fair_steep";
    (* Value of B *)
    b = bValueGen[regime][a];
    (* Range for y *)
    yMax = 4;
    (*
      ----------------------------------------------------------------
      Starting points along connected contour
        T = 1/2
      ----------------------------------------------------------------
    *)
    (* Value of T along contour *)
    tValue = 1/2;
    (* Seek starting point (x, y == 0) therealong *)
    yInit = 0;
    xInit = SeekRoot[
      tKnown[b][#, yInit] - tValue &,
      {0, yMax}
    ];
    (* (Probable) upper bound for arc length traversed *)
    sMax = 6;
    (* Contour (y > 0 half) *)
    xyContour =
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyContourSystem[b],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {
                tKnown[b][x[s], y[s]] < 0,
                x[s] < 0,
                y[s] > yMax
              },
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax},
          NoExtrapolation
        ]
      ];
    (* Actual arc length traversed *)
    sStart = DomainStart[xyContour];
    sEnd = DomainEnd[xyContour];
    (* Return starting points *)
    nMax = 12;
    startXYGen[a][regime]["connected"] = (
      Table[
        xyContour[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[sStart, sEnd, nMax]}]
        // Rest
        // Most
    );
    (*
      ----------------------------------------------------------------
      Hyperbolic critical terminal point (x = x_\[Sharp], y = 0)
      (The trivial point (x = x_\[Flat], y = 0) is useless.)
      ----------------------------------------------------------------
    *)
    startXYGen[a][regime]["hyperbolic"] = {
      {xSharp[a, b], 0}
    }
  ]
, {a, aValuesGen}];


(* ::Subsubsubsection:: *)
(*Steep regime B > 1*)


Table[
  Module[
   {yReflect,
    regime,
    b, yMax,
    tValue, xInit, yInit,
    sMax, xyContour,
    sStart, sEnd,
    nMax
   },
    (* Reflect y coordinate *)
    yReflect = # * {1, -1} &;
    (* Regime name *)
    regime = "steep";
    (* Value of B *)
    b = bValueGen[regime][a];
    (* Range for y *)
    yMax = 4;
    (*
      ----------------------------------------------------------------
      Starting points along connected contour
        T = 1/2
      ----------------------------------------------------------------
    *)
    (* Value of T along contour *)
    tValue = 1/2;
    (* Seek starting point (x, y == 0) therealong *)
    yInit = 0;
    xInit = SeekRoot[
      tKnown[b][#, yInit] - tValue &,
      {0, yMax}
    ];
    (* (Probable) upper bound for arc length traversed *)
    sMax = 6;
    (* Contour (y > 0 half) *)
    xyContour =
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyContourSystem[b],
            x[0] == xInit, y[0] == yInit,
            WhenEvent[
              {
                tKnown[b][x[s], y[s]] < 0,
                x[s] < 0,
                y[s] > yMax
              },
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax},
          NoExtrapolation
        ]
      ];
    (* Actual arc length traversed *)
    sStart = DomainStart[xyContour];
    sEnd = DomainEnd[xyContour];
    (* Return starting points *)
    nMax = 12;
    startXYGen[a][regime]["connected"] = (
      Table[
        xyContour[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[sStart, sEnd, nMax]}]
        // Rest
        // Most
    );
    (*
      ----------------------------------------------------------------
      Hyperbolic critical terminal point (x = x_\[Sharp], y = 0)
      (The point (x = x_\[Flat], y = 0) lies in the unphysical region T < 0.)
      ----------------------------------------------------------------
    *)
    startXYGen[a][regime]["hyperbolic"] = {
      {xSharp[a, b], 0}
    }
  ]
, {a, aValuesGen}];


(* ::Subsubsubsection:: *)
(*Starting point groups for each regime*)


idListGen =
  With[{d = "disconnected", c = "connected", h = "hyperbolic"},
    Association[
      "gentle" -> {d, c},
      "gentle_fair" -> {d, c, h},
      "fair" -> {d, c, h},
      "fair_steep" -> {c, h},
      "steep" -> {c, h}
    ]
  ];


(* ::Subsubsubsection:: *)
(*List of regimes*)


regimeListGen = Keys[idListGen];


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
(*By default the lower branch is used since this corresponds*)
(*to the candidate boundary for y > 0.*)


xTraDer[a_, b_, upperBranch_: False] := Function[{x, y},
  With[
   {p = p[b][x, y],
    q = q[b][x, y],
    f = f[a, b][x, y],
    vi = vi[a, b][x, y],
    sign = If[upperBranch, -1, +1]
   },
    Divide[
      p q + sign f Sqrt[vi],
      q^2 - f^2
    ]
  ] // Evaluate
] // Evaluate;


(* ::Subsubsection:: *)
(*Traced boundary curvature*)


(* ::Text:: *)
(*See (r5.20) (Page r5-5).*)


curTra[a_, b_, upperBranch_: False] := Function[{x, y},
  Module[{d, xDer, xDer2, cur},
    (* Abbreviation for y-derivative *)
    d = Dt[#, y, Constants -> {a, b}] &;
    (* x' and x'' *)
    (*
      Note: these are by default for the lower branch,
      since xTraDer uses the lower branch
     *)
    xDer = xTraDer[a, b, upperBranch][x, y];
    xDer2 = d[xDer];
    (* Curvature evaluated *)
    cur = xDer2 /. {d[x] -> xDer};
    cur = cur /. {d @ If[_, _, _] -> 0}
    (*
      Note: ugly d @ If[_, _, _] replacement needed
      since upperBranch returns an If expression.
      d[_If] doesn't work since it evaluates.
      d @ If[___] raises the warning If::argbu.
     *)
  ] // Evaluate
] // Evaluate;


(* ::Subsubsection:: *)
(*System of ODES for inflection frontier \[Kappa] = 0*)


curContourSystem[a_, b_, sSign_: 1] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{fun, p, q, slope},
      (* Scalar function \[Kappa] whose contours are sought *)
      (*
        Note: this is for the lower branch,
        since curTra uses the lower branch
       *)
      fun = curTra[a, b][x, y];
      (* Components of the gradient vector *)
      p = D[fun, x];
      q = D[fun, y];
      (* Magnitude of the gradient vector *)
      slope = Sqrt[p^2 + q^2];
      (* Return system of ODEs *)
      {
        x' == Re[q / slope],
        y' == Re[-p / slope]
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


(* ::Subsubsection:: *)
(*General case (B arbitrary)*)


(* ::Subsubsubsection:: *)
(*Candidate natural boundary (B = B_\[Natural](A), through x = x_\[Natural])*)


xTraCandNatGen[a_?NumericQ, terminateAtStraightContour_: False] :=
  With[{x = \[FormalX]},
    Module[{b, yMax},
      b = bNat[a];
      yMax = 5;
      NDSolveValue[
        {
          x'[y] == Re @ xTraDer[a, b][x[y], y],
          x[0] == xNat[a],
          WhenEvent[
            terminateAtStraightContour && x[y] > xStraight,
            "StopIntegration"
          ]
        }, x, {y, 0, yMax},
        NoExtrapolation
      ]
    ]
  ];


(* ::Subsection:: *)
(*Traced boundaries x = x(s), y = y(s)*)


(* ::Subsubsection:: *)
(*System of ODES for tracing*)


xyTraSystem[a_, b_, upperBranch_: True] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    With[
     {p = p[b][x, y],
      q = q[b][x, y],
      f = f[a, b][x, y],
      vi = vi[a, b][x, y],
      sign = If[upperBranch, +1, -1]
     },
      Module[{grad2, xDer, yDer},
        (* Square of gradient *)
        grad2 = p^2 + q^2;
        (* Return system of ODEs *)
        xDer = (-q f + sign p Re @ Sqrt[vi]) / grad2;
        yDer = (+p f + sign q Re @ Sqrt[vi]) / grad2;
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
              }, {x, y}, {s, -sMax, sMax},
              NoExtrapolation
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
xyTraCandSimp[a_?NumericQ, terminateAtStraightContour_: False] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{b, sMax},
      b = 1;
      sMax = 5;
      NDSolveValue[
        {
          xyTraSystem[a, b],
          x[0] == x0Simp[a], y[0] == 0,
          WhenEvent[
            terminateAtStraightContour && x[s] > xStraight,
            "StopIntegration"
          ]
        }, {x, y}, {s, 0, sMax},
        NoExtrapolation
      ]
    ]
  ];


(* ::Subsubsection:: *)
(*General case (B arbitrary)*)


(* ::Subsubsubsection:: *)
(*Generic boundaries (for the chosen representative values of A)*)


Module[
 {aValues,
  b, idList, sMax,
  xyInitList, xInit, yInit
 },
  (* Representative values of A *)
  aValues = aValuesGen;
  (* For each value of A *)
  Table[
    (* For each regime *)
    Table[
      (* Value of B *)
      b = bValueGen[regime][a];
      (* Starting point groups *)
      idList = idListGen[regime];
      (* Range for s *)
      sMax = 6;
      (* For each group of starting points *)
      Table[
        (* Starting points *)
        xyInitList = startXYGen[a][regime][id];
        (* Solve for traced boundaries *)
        xyTraGen[a][regime][id] =
          Table[
            {xInit, yInit} = xyInit;
            With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
              NDSolveValue[
                {
                  xyTraSystem[a, b],
                  x[0] == xInit, y[0] == yInit,
                  WhenEvent[
                    {
                      vi[a, b][x[s], y[s]] < 0,
                      tKnown[b][x[s], y[s]] < 0
                    },
                    "StopIntegration"
                  ]
                }, {x, y}, {s, -sMax, sMax},
                NoExtrapolation
              ]
            ]
          , {xyInit, xyInitList}]
      , {id, idList}]
    , {regime, regimeListGen}]
  , {a, aValues}]
];


(* ::Subsubsubsection:: *)
(*Constructing an asymmetric domain*)


(* ::Text:: *)
(*Non-trivial intersections (i.e. not (x_\[Flat], 0) or (x_\[Sharp], 0)) for the lower branch*)
(*  between the inflection frontier and the x-axis (y = 0), "axis inflection", and*)
(*  between the inflection frontier and the straight boundary (x = \[Pi]/2), "straight inflection".*)


Module[
 {a, b,
  xFl, xSh, xMid,
  sMax,
  xyInflAxis, xyInflStraight
 },
  a = aAsymm;
  b = bAsymm;
  (* Useful constants *)
  xFl = xFlat[a, b];
  xSh = xSharp[a, b];
  xMid = Way[xFl, xSh];
  (*
    ------------------------------------------------
    Axis inflection
    ------------------------------------------------
    xInflAxis "x_i(axis)":
      non-trivial x along y == 0 where inflection occurs
    yInflAxis "y_i(axis)":
      y where traced boundary through (x_i(axis), 0) intersects x == \[Pi]/2
   *)
  xInflAxis =
    SeekFirstRootBisection[
      curTra[a, b][#, 0] &,
      {xMid, xSh}
    ];
  sMax = 5;
  xyInflAxis =
    With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
      NDSolveValue[
        {
          xyTraSystem[a, b],
          x[0] == xInflAxis, y[0] == 0,
          WhenEvent[
            x[s] > xStraight,
            "StopIntegration"
          ]
        }, {x, y}, {s, 0, sMax},
        NoExtrapolation
      ]
    ];
  yInflAxis = -xyInflAxis[[2]] @ DomainEnd[xyInflAxis];
  (*
    ------------------------------------------------
    Straight inflection
    ------------------------------------------------
    yInflStraight "y_i(straight)":
      non-trivial y along x == \[Pi]/2 where inflection occurs
    xInflStraight "x_i(straight)":
      x where traced boundary through (\[Pi]/2, y_i(straight)) intersects y == 0
   *)
  yInflStraight =
    SeekFirstRootBisection[
      curTra[a, b][xStraight, #] &,
      {0.4, 0.8}
    ];
  sMax = 5;
  xyInflStraight =
    With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
      NDSolveValue[
        {
          xyTraSystem[a, b],
          x[0] == xStraight, y[0] == -yInflStraight,
          WhenEvent[
            y[s] > 0,
            "StopIntegration"
          ]
        }, {x, y}, {s, -sMax, 0},
        NoExtrapolation
      ]
    ];
  xInflStraight = xyInflStraight[[1]] @ DomainStart[xyInflStraight];
]


(* ::Text:: *)
(*The asymmetric domain shall have corner at*)
(*  x-coordinate 1/4 of the way from x_i(straight) to x_i(axis), and*)
(*  y-coordinate 9/10 of the way from 0 to the negative y, along that x,*)
(*    where inflection occurs for the lower branch.*)


Module[
 {a, b,
  yInflNegative,
  idList, upperBranch,
  sMax, sLower, sUpper
 },
  a = aAsymm;
  b = bAsymm;
  (*
    ------------------------------------------------
    Asymmetric domain corner
    ------------------------------------------------
    xAsymmCorner:
      1/4 of the way from x_i(straight) to x_i(axis)
    yAsymmCorner:
      9/10 of the way from 0 to the negative y, along that x,
      where inflection occurs for the lower branch.
   *)
  xAsymmCorner = Way[xInflStraight, xInflAxis, 1/4];
  yInflNegative = SeekFirstRootBisection[
    curTra[a, b][xAsymmCorner, #] &,
    {0, -0.5}
  ];
  yAsymmCorner = Way[0, yInflNegative, 9/10];
  (*
    ------------------------------------------------
    Asymmetric domain boundaries
    ------------------------------------------------
   *)
  idList = {"upper", "lower"};
  sMax = 2;
  Table[
    upperBranch = id === "upper";
    sLower = If[upperBranch, 0, -sMax];
    sUpper = If[upperBranch, sMax, 0];
    xyTraAsymm[id] =
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyTraSystem[a, b, upperBranch],
            x[0] == xAsymmCorner, y[0] == yAsymmCorner,
            WhenEvent[
              x[s] > xStraight,
              "StopIntegration"
            ]
          }, {x, y}, {s, sLower, sUpper},
          NoExtrapolation
        ]
      ];
  , {id, idList}];
];


(* ::Subsection:: *)
(*Numerical verification (finite elements)*)


(* ::Subsubsection:: *)
(*Simple case (B = 1)*)


(* ::Text:: *)
(*Here we consider the thin, lens-like convex domains.*)


(* ::Subsubsubsection:: *)
(*Generate finite element mesh*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
Module[
 {dest,
  yReflect,
  xyRad, xyRadEvenExt,
  xStart, sEnd, xEnd, yEnd,
  sSpacing, bPointList, nB,
  bMesh, mesh,
  prRad, prBath
 },
  Table[
    dest = FString[
      "cosine_simple-verification-mesh-{aNamesSimpConvex[a]}.txt"
    ];
    ExportIfNotExists[dest,
      (* Reflect y coordinate *)
      yReflect = # * {1, -1} &;
      (* Radiation boundary x == x(s), y == y(s) *)
      (* NOTE: y decreases with increasing s *)
      xyRad = xyTraCandSimp[a, True];
      (* Left-hand extremity thereof *)
      xStart = xyRad[[1]][0];
      (* Upper endpoint thereof *)
      sEnd = DomainEnd[xyRad];
      {xEnd, yEnd} = xyRad[sEnd] // Through // yReflect;
      (* Check that endpoint x value equals \[Pi]/2 *)
      If[xEnd == xStraight, Null,
        Print @ FString @ StringJoin[
          "WARNING (a == {N[a]}): ",
          "xEnd == {N[xEnd]} does not equal \[Pi]/2"
        ]
      ];
      (* Make spacing of boundary points 1/5 of the width *)
      sSpacing = 1/5 (xEnd - xStart);
      (* Boundary points *)
      xyRadEvenExt[s_] /; s >= 0 := xyRad[s] // Through;
      xyRadEvenExt[s_] /; s < 0 := xyRadEvenExt[-s] // yReflect;
      bPointList = Join[
        (* Radiation boundary (x(s), y(s)) *)
        Table[
          xyRadEvenExt[s]
        , {s, UniformRange[-sEnd, sEnd, sSpacing]}],
        (* Straight boundary x == \[Pi]/2 *)
        Table[
          {xStraight, y}
        , {y, UniformRange[-yEnd, yEnd, sSpacing]}]
          // Rest // Most
      ];
      bPointList = bPointList;
      nB = Length[bPointList];
      (* Build boundary element mesh *)
      bMesh = ToBoundaryMesh[
        "Coordinates" -> bPointList,
        "BoundaryElements" -> {
          LineElement[
            Table[{n, n + 1}, {n, nB}] // Mod[#, nB, 1] &
          ]
        }
      ];
      (* Build mesh *)
      mesh = ToElementMesh[bMesh,
        "ImproveBoundaryPosition" -> True
      ];
      (* Predicate functions for radiation and bath (straight) boundaries *)
      prRad = Function[{x, y}, x < xStraight // Evaluate];
      prBath = Function[{x, y}, x == xStraight // Evaluate];
      (* Export *)
      {a, mesh, prRad, prBath} // Compress
    ]
  , {a, aValuesSimpConvex}]
]


(* ::Subsubsubsection:: *)
(*Solve boundary value problem (Version 12 required)*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
Table[
  Module[
   {dest, source,
    a, mesh, prRad, prBath,
    tSol
   },
    dest = FString[
      "cosine_simple-verification-solution-{aNamesSimpConvex[a]}.txt"
    ];
    ExportIfNotExists[dest,
      (* Import mesh *)
      source = FString[
        "cosine_simple-verification-mesh-{aNamesSimpConvex[a]}.txt"
      ];
      {a, mesh, prRad, prBath} = Import[source] // Uncompress;
      (* Solve boundary value problem *)
      With[{x = \[FormalX], y = \[FormalY], t = \[FormalCapitalT]},
        tSol = NDSolveValue[
          {
            (* Steady state heat equation *)
            -Laplacian[t[x, y], {x, y}] ==
              (* Radiation boundary condition *)
              NeumannValue[(-1 / a) t[x, y]^4, prRad[x, y]],
            (* Heat bath (straight edge) boundary condition *)
            DirichletCondition[t[x, y] == 1, prBath[x, y]]
          }, t, Element[{x, y}, mesh]
        ]
      ];
      tSol // Compress
    ]
  ]
, {a, aValuesSimpConvex}]


(* ::Subsubsection:: *)
(*General case (B arbitrary)*)


(* ::Text:: *)
(*Here we consider the asymmetric convex domain constructed above.*)


(* ::Subsubsubsection:: *)
(*Generate finite element mesh*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
Module[
 {dest,
  xyUpper, xyLower, yTop, yBottom,
  sSpacing, bPointList, nB,
  bMesh, mesh,
  prRad, prBath
 },
  dest = "cosine_general-verification-mesh-asymmetric.txt";
  ExportIfNotExists[dest,
    (* Top- and bottom-corner y-coordinates *)
    (*
      NOTE: "upper" and "lower" refer to the branch,
      not the relative vertical position of the two boundary curves
     *)
    xyUpper = xyTraAsymm["upper"];
    xyLower = xyTraAsymm["lower"];
    yTop = xyLower[[2]] @ DomainStart[xyLower];
    yBottom = xyUpper[[2]] @ DomainEnd[xyUpper];
    (* Make spacing of boundary points 1/50 of the height *)
    sSpacing = 1/50 (yTop - yBottom);
    (* Boundary points *)
    bPointList = Join[
      (* Radiation boundary: lower branch (physically higher) *)
      Table[
        xyLower[s] // Through
      , {s, UniformRange[
          DomainStart[xyLower],
          DomainEnd[xyLower],
          sSpacing
        ]
      }] // Most,
      (* Radiation boundary: upper branch (physically lower) *)
      Table[
        xyUpper[s] // Through
      , {s, UniformRange[
          DomainStart[xyUpper],
          DomainEnd[xyUpper],
          sSpacing
        ]
      }] // Most,
      (* Straight boundary x == \[Pi]/2 *)
      Table[
        {xStraight, y}
      , {y, UniformRange[yBottom, yTop, sSpacing]}] // Most
    ];
    nB = Length[bPointList];
    (* Build boundary element mesh *)
    bMesh = ToBoundaryMesh[
      "Coordinates" -> bPointList,
      "BoundaryElements" -> {
        LineElement[
          Table[{n, n + 1}, {n, nB}] // Mod[#, nB, 1] &
        ]
      }
    ];
    (* Build mesh *)
    mesh = ToElementMesh[bMesh,
      "ImproveBoundaryPosition" -> True
    ];
    (* Predicate functions for radiation and bath (straight) boundaries *)
    prRad = Function[{x, y}, x < xStraight // Evaluate];
    prBath = Function[{x, y}, x == xStraight // Evaluate];
    (* Export *)
    {mesh, prRad, prBath} // Compress
  ]
]


(* ::Subsubsubsection:: *)
(*Solve boundary value problem (Version 12 required)*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
Module[
 {dest, source,
  mesh, prRad, prBath,
  a, tSol
 },
  dest = "cosine_general-verification-solution-asymmetric.txt";
  ExportIfNotExists[dest,
    (* Import mesh *)
    source = "cosine_general-verification-mesh-asymmetric.txt";
    {mesh, prRad, prBath} = Import[source] // Uncompress;
    (* Solve boundary value problem *)
    a = aAsymm;
    With[{x = \[FormalX], y = \[FormalY], t = \[FormalCapitalT]},
      tSol = NDSolveValue[
        {
          (* Steady state heat equation *)
          -Laplacian[t[x, y], {x, y}] ==
            (* Radiation boundary condition *)
            NeumannValue[(-1 / a) t[x, y]^4, prRad[x, y]],
          (* Heat bath (straight edge) boundary condition *)
          DirichletCondition[t[x, y] == 1, prBath[x, y]]
        }, t, Element[{x, y}, mesh]
      ]
    ];
    tSol // Compress
  ]
]


(* ::Subsection:: *)
(*Italicise symbols*)


aIt = Italicise["A"];
bIt = Italicise["B"];
sIt = Italicise["s"];
xIt = Italicise["x"];
yIt = Italicise["y"];


(* ::Subsection:: *)
(*Options for exported GIFS*)


gifOpts = Sequence[
  AnimationRepetitions -> Infinity,
  "DisplayDurations" -> 0.5
];


(* ::Subsection:: *)
(*Global styles for plots*)


textStyle = Style[#, 18] &;


natStyle = Red;
flatStyle = Blue;
sharpStyle = Orange;


contStyle = LightGray;
straightStyle = Directive[Dashed, Magenta];
nonStyle = Directive[Opacity[0.7], LightGray];
termStyle = Pink;
unphysStyle = Black;
inflStyle = Directive[Thick, Green];


upperStyle = Blue;
lowerStyle = Red;
convexStyle = Black;


simp0Style = Blue;


pointStyle = PointSize[Large];
startingPointStyle = Directive[Opacity[0.7], pointStyle];
glowStyle = Directive[Thick, Yellow, Opacity[0.7]];


inflDotStyle = Directive[Red, Opacity[0.7], pointStyle];


(* ::Subsection:: *)
(*Computations for figures*)


(* ::Subsubsection:: *)
(*Simple case (B = 1) traced boundaries*)


(* ::Subsubsubsection:: *)
(*Starting points along T = 0 contour*)


startXYSimpFigure =
  Module[{b, yMax, num, eps, yValues},
    (* Constants *)
    b = 1;
    yMax = 2;
    num = 16;
    eps = 10^-4;
    yValues = Subdivide[-yMax, yMax, num];
    (* Return starting points *)
    Table[
      {ArcCos[(1 - eps) Sech[y] / b], y}
    , {y, yValues}]
  ];


(* ::Subsubsubsection:: *)
(*Traced boundaries*)


xyTraSimpFigure =
  Module[{a, b, sMax, viTol, xInit, yInit},
    a = 1/2;
    b = 1;
    sMax = 4;
    viTol = 10^-6;
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
    , {xyInit, startXYSimpFigure}]
  ];


(* ::Subsubsection:: *)
(*Simple case (B = 1) traced boundaries, patched spiky*)


(*
  Corners (x_1, y_1), ..., (x_n, y_n) are (carefully) chosen,
  then the intersections between lower-branch(i) and upper-branch(i+1),
  for i == 1, ..., n - 1, chosen. Note that within each spike,
  the lower-branch curve is physically higher (in y).
*)


patchedCornerList = {
  {1, -0.5},
  {0.7, -0.25},
  {0.85, -0.13},
  {1.1, 0.2},
  {1.05, 0.3},
  {0.9, 0.8},
  {1.3, 1.23},
  Nothing
};


patchedCornerNum = Length[patchedCornerList];


patchedBoundaryUpperList =
  Module[{a, b, sMax, viTol, xInit, yInit},
    a = 1/2;
    b = 1;
    sMax = 4;
    viTol = 0;
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
    , {xyInit, patchedCornerList}]
  ];


patchedBoundaryLowerList =
  Module[{a, b, sMax, viTol, xInit, yInit},
    a = 1/2;
    b = 1;
    sMax = 4;
    viTol = 0;
    Table[
      {xInit, yInit} = xyInit;
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            xyTraSystem[a, b, False],
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
    , {xyInit, patchedCornerList}]
  ];


patchedIntersectionList =
  Table[
    SeekParametricIntersection[
      patchedBoundaryLowerList[[i]],
      patchedBoundaryUpperList[[i + 1]]
    ]
  , {i, patchedCornerNum - 1}];


{patchedIntersectionLowerList, patchedIntersectionUpperList} =
  Transpose @ patchedIntersectionList;
AppendTo[patchedIntersectionLowerList,
  DomainStart @ patchedBoundaryLowerList[[patchedCornerNum]]
];
PrependTo[patchedIntersectionUpperList,
  DomainEnd @ patchedBoundaryUpperList[[1]]
];


(* ::Subsubsection:: *)
(*Simple case (B = 1) traced boundaries, patched smooth*)


patchedBoundarySmooth =
  Module[{a, b, sMax, viTol},
    a = 1/2;
    b = 1;
    sMax = 4;
    viTol = 0;
    With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
      NDSolveValue[
        {
          xyTraSystem[a, b],
          x[0] == x0Simp[a], y[0] == 0,
          WhenEvent[
            tKnown[b][x[s], y[s]] < 0,
            "StopIntegration"
          ]
        }, {x, y}, {s, -sMax, sMax},
        NoExtrapolation
      ]
    ]
  ];


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
  , {{b, bInit, bIt}, bMin, bMax, Appearance -> "Open"}]
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
        ContourLabels -> None,
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
  , {{b, bInit, bIt}, bMin, bMax, Appearance -> "Open"}]
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
  numTo1, numBeyond1,
  resetLabel,
  critTermPoint
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
  (* Label for reset *)
  resetLabel[var_] := Row @ {"Reset ", var};
  (* Critical terminal point *)
  critTermPoint[x0_, style_] :=
    Graphics @ {Directive[style, pointStyle],
      Point @ {x0, 0}
    };
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
        ContourLabels -> None,
        Contours -> numTo1 + numBeyond1,
        ContourShading -> None,
        ContourStyle -> contStyle,
        PlotRange -> {0, 1 + (1 + numBeyond1) / numTo1}
      ],
      (* Straight contour *)
      Graphics @ {straightStyle,
        Line @ {{xStraight, -yMaxCont}, {xStraight, yMaxCont}}
      },
      (* Critical terminal points along y == 0 *)
      Which[
        (* Two distinct terminal points, x_\[Flat] & x_\[Sharp] *)
        b > bNat[a],
        {
          critTermPoint[Re @ xFlat[a, b], flatStyle],
          critTermPoint[Re @ xSharp[a, b], sharpStyle]
          (* Re chops off small imaginary part in floating point arithmetic *)
        },
        (* One terminal point x_\[Natural] *)
        b == bNat[a],
        {
          critTermPoint[xNat[a], natStyle]
        },
        (* Zero critical terminal points *)
        True,
        {}
      ]
    ]
  , {{a, aInit, aIt}, aMin, aMax, Appearance -> "Open"}
  , {{b, bInit, bIt}, bMin, bMax, Appearance -> "Open"}
  , Button[resetLabel[aIt], a = aInit]
  , Button[resetLabel[bIt], b = bInit]
  ]
]


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
  , {{a, aInit, aIt}, aMin, aMax, Appearance -> "Open"}
  , {{b, bInit, bIt}, bMin, bMax, Appearance -> "Open"}]
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
      ContourLabels -> None,
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
        ContourLabels -> None,
        ContourStyle -> simp0Style,
        PlotRange -> {0, 1}
      ]
    ]
  , {a, aValues}]
] // Ex["cosine_simple-viable.gif", gifOpts]


(* ::Subsection:: *)
(*Animation for viable domain, general case (A = const, B varying)*)


Module[
 {a,
  bStep, bCrit, bMin, bMax, bValues,
  xMin, xMax, yMax,
  mar,
  xMinMar, xMaxMar, yMaxMar,
  straightContour,
  critTermPoint
 },
  (* Value of A *)
  a = 3;
  (* Values of B *)
  bStep = 1/100;
  bCrit = bNat[a];
  bMin = Max[0, Floor[bCrit - 4 bStep, bStep]];
  bMax = Ceiling[bCrit + 8 bStep, bStep];
  bValues = Range[bMin, bMax, bStep];
  bValues = Append[bCrit] @ bValues;
  bValues = bValues // Sort[#, Less] &;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 5/4;
  yMax = 2;
  (* Plot ranges with margin *)
  mar = 1/5 xMax;
  xMinMar = xMin - mar;
  xMaxMar = xMax + mar;
  yMaxMar = yMax + mar;
  (* Straight contour (independent of A and B) *)
  straightContour =
    Graphics @ {straightStyle,
      Line @ {{xStraight, -yMaxMar}, {xStraight, yMaxMar}}
    };
  (* Critical terminal point *)
  critTermPoint[x0_, style_] :=
    Graphics @ {Directive[style, pointStyle],
      Point @ {x0, 0}
    };
  (* Animation *)
  Table[
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        ImageSize -> 180,
        PlotLabel -> BoxedLabel[
          Row[
            {
              aIt == a,
              bIt == N[b, 2]
            } // If[b == bCrit,
              Style[#, natStyle] & /@ # &,
              Identity
            ],
            ","
          ],
          FrameStyle -> If[b == bCrit,
            natStyle,
            Automatic
          ]
        ]
      ],
      (* Unphysical domain *)
      RegionPlot[tKnown[b][x, y] < 0,
        {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
        BoundaryStyle -> unphysStyle,
        PlotPoints -> 50,
        PlotStyle -> unphysStyle
      ],
      (* Non-viable domain *)
      RegionPlot[vi[a, b][x, y] < 0,
        {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
        BoundaryStyle -> termStyle,
        PlotPoints -> 50,
        PlotStyle -> nonStyle
      ],
      (* Known solution contours *)
      ContourPlot[tKnown[b][x, y],
        {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
        ContourLabels -> None,
        ContourShading -> None,
        ContourStyle -> contStyle,
        PlotRange -> {0, tKnown[b][xMax, 0]}
      ],
      straightContour,
      (* Critical terminal points along y == 0 *)
      Which[
        (* Two distinct terminal points, x_\[Flat] & x_\[Sharp] *)
        b > bCrit,
        {
          critTermPoint[Re @ xFlat[a, b], flatStyle],
          critTermPoint[Re @ xSharp[a, b], sharpStyle]
          (* Re chops off small imaginary part in floating point arithmetic *)
        },
        (* One terminal point x_\[Natural] *)
        b == bCrit,
        {
          critTermPoint[N @ xNat[a], natStyle]
          (* Doesn't work without applying N for some reason *)
        },
        (* Zero critical terminal points *)
        True,
        {}
      ]
    ]
  , {b, bValues}]
] // Ex["cosine_general-viable.gif", gifOpts]


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
      ContourLabels -> None,
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
    xCand, yCand
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
          {xCand, yCand} = xyTraCandSimp[a];
          (* Evaluate residual therealong *)
          res[type] /. {
            x' -> xCand'[s],
            y' -> yCand'[s],
            x -> xCand[s],
            y -> yCand[s]
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


(* ::Text:: *)
(*More compact version:*)


With[{y = \[FormalY], a = \[FormalCapitalA], c = \[FormalCapitalC], s = \[FormalCapitalS]},
  Block[{$Assumptions = 0 < a <= 1 && 0 < y < ArcCosh[1/a]},
    curTra[a, 1][Pi/2, y] == (
      a^2 c / Sqrt[a^2 c^2 - 1] * (
        2 s - (1 + s^2) (a^2 s + 4 Sqrt[a^2 (1 + s^2) - 1])
      )
    )
      /. {c -> Cosh[y], s -> Sinh[y]}
      // FullSimplify
  ]
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


(* ::Section:: *)
(*General case (B arbitrary) plots*)


(* ::Subsection:: *)
(*Critical terminal points along y = 0*)


Module[
 {a,
  bCrit, xCrit,
  bMin, bMax,
  xMin, xMax
 },
  (* Value of A *)
  a = 3;
  (* Critical values B_\[Natural] and x_\[Natural] *)
  (* (N needed otherwise Epilog doesn't work) *)
  bCrit = bNat[a] // N;
  xCrit = xNat[a] // N;
  (* B range *)
  bMin = 0;
  bMax = 5;
  (* x range *)
  xMin = 0;
  xMax = 1.1 xSharp[a, #] & [
    With[{b = \[FormalCapitalB]},
      b /. Last @ Solve[
        D[xSharp[a, b], b] == 0,
        b
      ]
    ]
  ];
  (* Plot *)
  Plot[
    (* x_\[Sharp] & x_\[Flat]*)
    {
      xSharp[a, b],
      xFlat[a, b],
      Piecewise @ {
        {xMin, b < bCrit},
        {Indeterminate, True}
      }
    },
    {b, bMin + 10^-6, bMax},
    AxesLabel -> {bIt, xIt},
    Epilog -> {
      (* x == x_\[Natural], B == B_\[Natural] *)
      Directive[pointStyle, natStyle],
      Point @ {bCrit, xCrit},
      Text[
        Column[
          textStyle /@ {
            Subscript[xIt, "nat"] == xCrit,
            Subscript[bIt, "nat"] == bCrit
          }
        ],
        {bCrit, xCrit},
        {-1.2, 0}
      ],
    },
    Filling -> {
      1 -> Top,
      2 -> Bottom,
      3 -> Top
    },
    FillingStyle -> nonStyle,
    PlotLabel -> BoxedLabel[aIt == a],
    PlotLegends -> (
      xIt == Subscript[xIt, #] & /@
        {"sharp", "flat"}
    ),
    PlotRange -> {xMin, xMax},
    PlotStyle -> {flatStyle, sharpStyle, None},
    PlotOptions[Axes] // Evaluate
  ] // PrettyString[
    "flat" -> "\[Flat]",
    "nat" -> "\[Natural]",
    "sharp" -> "\[Sharp]"
  ] // Ex @ FString[
    "cosine_general-critical-a-{ToName[a]}.pdf"
  ]
]


(* ::Subsection:: *)
(*Starting points for boundary tracing*)


(* (This is somewhat slow, as 10 files are being generated in one go.) *)
Module[
 {xMin, xMax, yMax,
  tContourNum, tContourRange,
  mar,
  xMinMar, xMaxMar, yMaxMar,
  emptyFrame,
  unphysicalDomain,
  nonViableDomain,
  generalContours,
  straightContour,
  labelFun,
  regime, b, idList
 },
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 3/2;
  yMax = 4;
  tContourNum = 13;
  tContourRange = {0, 2};
  (* Plot range with margin *)
  mar = 1/5 xMax;
  xMinMar = xMin - mar;
  xMaxMar = xMax + mar;
  yMaxMar = yMax + mar;
  (* Empty frame *)
  emptyFrame[a_, b_] := (
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 240,
      PlotLabel -> BoxedLabel[
        Row[
          {aIt == a, bIt == N[b]},
          ","
        ]
      ]
    ]
  );
  (* Unphysical domain *)
  unphysicalDomain[b_] := (
    RegionPlot[tKnown[b][x, y] < 0,
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      BoundaryStyle -> unphysStyle,
      PlotPoints -> 50,
      PlotStyle -> unphysStyle
    ]
  );
  (* Non-viable domain *)
  nonViableDomain[a_, b_] := (
    RegionPlot[vi[a, b][x, y] < 0,
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      BoundaryStyle -> termStyle,
      PlotPoints -> 70,
      PlotStyle -> nonStyle
    ]
  );
  (* Known solution contours (general) *)
  generalContours[b_] := (
    ContourPlot[
      tKnown[b][x, y],
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      ContourLabels -> None,
      Contours -> tContourNum,
      ContourShading -> None,
      ContourStyle -> contStyle,
      PlotRange -> tContourRange
    ]
  );
  (* Known solution contours (straight) *)
  straightContour = Graphics @ {straightStyle,
    Line @ {{xStraight, -yMaxMar}, {xStraight, yMaxMar}}
  };
  (* Labelling function *)
  labelFun = Function @ Placed[
    #2[[2]],
    Center
  ];
  (* For each value of A *)
  Table[
    (* For each regime *)
    Table[
      (* Value of B *)
      b = bValueGen[regime][a];
      (* Starting point groups *)
      idList = idListGen[regime];
      (* Plot *)
      Show[
        emptyFrame[a, b],
        unphysicalDomain[b],
        nonViableDomain[a, b],
        generalContours[b],
        straightContour,
        ListPlot[
          Table[
            startXYGen[a][regime][id]
          , {id, idList}],
          LabelingFunction -> labelFun,
          PlotLegends -> idList,
          PlotStyle -> startingPointStyle
        ]
      ] // Ex @ FString[
        "cosine_general-a-{ToName[a]}-{regime}-traced-starting.pdf"
      ]
    , {regime, regimeListGen}]
  , {a, aValuesGen}]
]


(* ::Subsection:: *)
(*Traced boundaries*)


(* (This is somewhat slow, as 10 files are being generated in one go.) *)
Module[
 {xMin, xMax, yMax,
  tContourNum, tContourRange,
  mar,
  xMinMar, xMaxMar, yMaxMar,
  emptyFrame,
  unphysicalDomain,
  nonViableDomain,
  generalContours,
  straightContour,
  regime, b, idList
 },
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 3/2;
  yMax = 4;
  tContourNum = 13;
  tContourRange = {0, 2};
  (* Plot range with margin *)
  mar = 1/5 xMax;
  xMinMar = xMin - mar;
  xMaxMar = xMax + mar;
  yMaxMar = yMax + mar;
  (* Empty frame *)
  emptyFrame[a_, b_] := (
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 240,
      PlotLabel -> BoxedLabel[
        Row[
          {aIt == a, bIt == N[b]},
          ","
        ]
      ]
    ]
  );
  (* Unphysical domain *)
  unphysicalDomain[b_] := (
    RegionPlot[tKnown[b][x, y] < 0,
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      BoundaryStyle -> unphysStyle,
      PlotPoints -> 50,
      PlotStyle -> unphysStyle
    ]
  );
  (* Non-viable domain *)
  nonViableDomain[a_, b_] := (
    RegionPlot[vi[a, b][x, y] < 0,
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      BoundaryStyle -> termStyle,
      PlotPoints -> 70,
      PlotStyle -> nonStyle
    ]
  );
  (* Known solution contours (general) *)
  generalContours[b_] := (
    ContourPlot[
      tKnown[b][x, y],
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      ContourLabels -> None,
      Contours -> tContourNum,
      ContourShading -> None,
      ContourStyle -> contStyle,
      PlotRange -> tContourRange
    ]
  );
  (* Known solution contours (straight) *)
  straightContour = Graphics @ {straightStyle,
    Line @ {{xStraight, -yMaxMar}, {xStraight, yMaxMar}}
  };
  (* For each value of A *)
  Table[
    (* For each regime *)
    Table[
      (* Value of B *)
      b = bValueGen[regime][a];
      (* Starting point groups *)
      idList = idListGen[regime];
      (* Plot *)
      Show[
        emptyFrame[a, b],
        unphysicalDomain[b],
        nonViableDomain[a, b],
        generalContours[b],
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
          , {xy, xyTraGen[a][regime][id]}]
        , {id, idList}]
        (* Not bothering with hyperbolic traced boundaries *)
      ] // Ex @ FString[
        "cosine_general-a-{ToName[a]}-{regime}-traced-full.pdf"
      ]
    , {regime, regimeListGen}]
  , {a, aValuesGen}]
]


(* ::Subsection:: *)
(*Traced boundaries from x-axis*)


(* ::Subsubsection:: *)
(*Visualiser (poor man's Manipulate)*)


Module[
 {a, b,
  xMin, xMax, yMax,
  tContourNum, tContourRange,
  mar,
  xMinMar, xMaxMar, yMaxMar,
  emptyFrame,
  unphysicalDomain,
  nonViableDomain,
  generalContours,
  straightContour,
  lineOfSymmetry,
  sMax,
  xInitMin, xInitMax, xInitList,
  xyList,
  tracedBoundaries,
  xyMid, sInflList,
  xInitInfl, yInitInfl,
  xyInflList,
  inflectionFrontiers
 },
  (* Values of A and B (to be set manually) *)
  a = 12;
  b = 1.05 bNat[a];
  (*
    ------------------------------------------------
    Barely an infinite family of convex domains
    ------------------------------------------------
    a = 7.5;
    b = 1.05 bNat[a];
    ------------------------------------------------
    Probably only one convex domain (lens-like)
    ------------------------------------------------
    a = 7;
    b = 1.1 bNat[a];
    ------------------------------------------------
    Weird inflection frontier
    ------------------------------------------------
    a = 7.5;
    b = 1.3 bNat[a];
    ------------------------------------------------
    Slightly asymmetric convex domains
    ------------------------------------------------
    a = 12;
    b = 1.05 bNat[a];
   *)
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 3/2;
  yMax = 2;
  tContourNum = 13;
  tContourRange = {0, 2};
  (* Plot range with margin *)
  mar = 1/5 xMax;
  xMinMar = xMin - mar;
  xMaxMar = xMax + mar;
  yMaxMar = yMax + mar;
  (* Empty frame *)
  emptyFrame[a_, b_] := (
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 360,
      PlotLabel -> BoxedLabel[
        Row[
          {aIt == a, bIt == N[b]},
          ","
        ]
      ]
    ]
  );
  (* Unphysical domain *)
  unphysicalDomain[b_] := (
    RegionPlot[tKnown[b][x, y] < 0,
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      BoundaryStyle -> unphysStyle,
      PlotPoints -> 50,
      PlotStyle -> unphysStyle
    ]
  );
  (* Non-viable domain *)
  nonViableDomain[a_, b_] := (
    RegionPlot[vi[a, b][x, y] < 0,
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      BoundaryStyle -> termStyle,
      PlotPoints -> 70,
      PlotStyle -> nonStyle
    ]
  );
  (* Known solution contours (general) *)
  generalContours[b_] := (
    ContourPlot[
      tKnown[b][x, y],
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      ContourLabels -> None,
      Contours -> tContourNum,
      ContourShading -> None,
      ContourStyle -> contStyle,
      PlotRange -> tContourRange
    ]
  );
  (* Known solution contours (straight) *)
  straightContour = Graphics @ {straightStyle,
    Line @ {{xStraight, -yMaxMar}, {xStraight, yMaxMar}}
  };
  (* Line of symmetry (y == 0) *)
  lineOfSymmetry[a_, b_] := Graphics @ {Directive[Thin, Orange],
    Line @ N @ {{xFlat[a, b], 0}, {xSharp[a, b], 0}}
  };
  (* Range for s *)
  sMax = 6;
  (* Initial x values along axis *)
  xInitMin = xFlat[a, b];
  xInitMax = xSharp[a, b];
  xInitList = Subdivide[xInitMin, xInitMax, 4];
  (* Exclude x_\[Flat] == 0 for B = 1 *)
  If[b == 1, xInitList = xInitList // Rest];
  (* Compute traced boundaries *)
  xyList = Table[
    With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
      NDSolveValue[
        {
          xyTraSystem[a, b],
          x[0] == xInit, y[0] == 0,
          WhenEvent[
            {
              vi[a, b][x[s], y[s]] < 0,
              tKnown[b][x[s], y[s]] < 0
            },
            "StopIntegration"
          ]
        }, {x, y}, {s, -sMax, sMax},
        NoExtrapolation
      ]
    ]
  , {xInit, xInitList}];
  tracedBoundaries =
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // {#, {#[[1]], -#[[2]]}} &
          // Evaluate,
        {s, DomainStart[xy], DomainEnd[xy]},
        PlotStyle -> {upperStyle, lowerStyle}
      ]
    , {xy, xyList}];
  (* Compute inflection frontiers for the lower branch *)
  xyMid = xyList[[3]];
  sInflList =
    Table[
      SeekFirstRootBisection[
        curTra[a, b] @@ Through @ xyMid[#] &,
        {0, end[xyMid]}
      ]
    , {end, {DomainStart, DomainEnd}}];
  xyInflList =
    Table[
      {xInitInfl, yInitInfl} = xyMid[sInit] // Through;
      With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
        NDSolveValue[
          {
            curContourSystem[a, b],
            x[0] == xInitInfl, y[0] == yInitInfl,
            WhenEvent[
              {
                x[s] < xMinMar,
                x[s] > xMaxMar,
                Abs @ y[s] > yMaxMar,
                vi[a, b][x[s], y[s]] < 0
              },
              "StopIntegration"
            ]
          }, {x, y}, {s, -sMax, sMax},
          NoExtrapolation
        ]
      ]
    , {sInit, sInflList}];
  inflectionFrontiers[
    opt : OptionsPattern @ {"Mirror" -> False}
  ] :=
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // If[OptionValue["Mirror"],
            {#, {#[[1]], -#[[2]]}} &,
            Identity
          ]
          // Evaluate,
        {s, DomainStart[xy], DomainEnd[xy]},
        PlotStyle -> inflStyle
      ]
    , {xy, xyInflList}];
  (* Plot *)
  Show[
    emptyFrame[a, b],
    unphysicalDomain[b],
    nonViableDomain[a, b],
    generalContours[b],
    straightContour,
    lineOfSymmetry[a, b],
    tracedBoundaries,
    inflectionFrontiers[]
  ]
]


(* ::Subsubsection:: *)
(*Animation*)


Module[
 {bMin, bMax, bValues,
  xMin, xMax, yMax,
  tContourNum, tContourRange,
  mar,
  xMinMar, xMaxMar, yMaxMar,
  emptyFrame,
  unphysicalDomain,
  nonViableDomain,
  generalContours,
  straightContour,
  sMax,
  xInitMin, xInitMax, xInitList,
  xyList
 },
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 3/2;
  yMax = 2.5;
  tContourNum = 13;
  tContourRange = {0, 2};
  (* Plot range with margin *)
  mar = 1/5 xMax;
  xMinMar = xMin - mar;
  xMaxMar = xMax + mar;
  yMaxMar = yMax + mar;
  (* Empty frame *)
  emptyFrame[a_, b_] := (
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 200,
      PlotLabel -> BoxedLabel[
        Column[
          {aIt == a, bIt == N[b, 2]},
          Alignment -> Center
        ]
      ]
    ]
  );
  (* Unphysical domain *)
  unphysicalDomain[b_] := (
    RegionPlot[tKnown[b][x, y] < 0,
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      BoundaryStyle -> unphysStyle,
      PlotPoints -> 50,
      PlotStyle -> unphysStyle
    ]
  );
  (* Non-viable domain *)
  nonViableDomain[a_, b_] := (
    RegionPlot[vi[a, b][x, y] < 0,
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      BoundaryStyle -> termStyle,
      PlotPoints -> 70,
      PlotStyle -> nonStyle
    ]
  );
  (* Known solution contours (general) *)
  generalContours[b_] := (
    ContourPlot[
      tKnown[b][x, y],
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      ContourLabels -> None,
      Contours -> tContourNum,
      ContourShading -> None,
      ContourStyle -> contStyle,
      PlotRange -> tContourRange
    ]
  );
  (* Known solution contours (straight) *)
  straightContour = Graphics @ {straightStyle,
    Line @ {{xStraight, -yMaxMar}, {xStraight, yMaxMar}}
  };
  (* Range for s *)
  sMax = 6;
  (* For various values of A *)
  Table[
    (* Values of B *)
    bMin = bNat[a];
    bMax = 1;
    bValues = Subdivide[bMin, bMax, 20];
    (* B-manipulation *)
    Table[
      (* Initial x values along axis *)
      xInitMin = xFlat[a, b];
      xInitMax = xSharp[a, b];
      xInitList = Subdivide[xInitMin, xInitMax, 4];
      (* Exclude x_\[Flat] == 0 for B = 1 *)
      If[b == 1, xInitList = xInitList // Rest];
      (* Compute traced boundaries *)
      xyList = Table[
        With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
          NDSolveValue[
            {
              xyTraSystem[a, b],
              x[0] == xInit, y[0] == 0,
              WhenEvent[
                {
                  vi[a, b][x[s], y[s]] < 0,
                  tKnown[b][x[s], y[s]] < 0
                },
                "StopIntegration"
              ]
            }, {x, y}, {s, -sMax, sMax},
            NoExtrapolation
          ]
        ]
      , {xInit, xInitList}];
      (* Plot *)
      Show[
        emptyFrame[a, b],
        unphysicalDomain[b],
        nonViableDomain[a, b],
        generalContours[b],
        straightContour,
        (* Traced boundaries *)
        Table[
          ParametricPlot[
            xy[s]
              // Through
              // {#, {#[[1]], -#[[2]]}} &
              // Evaluate,
            {s, DomainStart[xy], DomainEnd[xy]},
            PlotStyle -> {upperStyle, lowerStyle}
          ]
        , {xy, xyList}]
      ]
    , {b, bValues}]
      // Ex[
        FString["cosine_general-a-{a}-traced-from-axis.gif"],
        gifOpts
      ]
  , {a, aValuesGen}]
]


(* ::Subsection:: *)
(*Curvature along x-axis*)


(* ::Subsubsection:: *)
(*Interactive visualiser*)


DynamicModule[
 {aInit, aMin, aMax,
  bInit, bMin, bMax,
  xMin, xMax,
  curMax, curMaxMar,
  resetLabel,
  critTermPoint
 },
  (* Values of A *)
  aInit = 12;
  aMin = 1;
  aMax = 20;
  (* Values of B *)
  bInit = 1.05 bNat[aInit];
  bMin = 0;
  bMax = 1;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 5/4;
  curMax = 1;
  curMaxMar = 1.2 curMax;
  (* Label for reset *)
  resetLabel[var_] := Row @ {"Reset ", var};
  (* Critical terminal point *)
  critTermPoint[x0_, style_] :=
    Graphics @ {Directive[style, pointStyle],
      Point @ {x0, 0}
    };
  (* Plot *)
  Manipulate[
    Show[
      Plot[
        curTra[a, b][x, 0], {x, xMin, xMax},
        AxesLabel -> {Automatic, "\[Kappa]"},
        PlotRange -> {-curMax, curMax},
        PlotOptions[Axes] // Evaluate
      ],
      (* Straight contour *)
      Graphics @ {straightStyle,
        Line @ {{xStraight, -curMaxMar}, {xStraight, curMaxMar}}
      },
      (* Critical terminal points along y == 0 *)
      Which[
        (* Two distinct terminal points, x_\[Flat] & x_\[Sharp] *)
        b > bNat[a],
        {
          critTermPoint[Re @ xFlat[a, b], flatStyle],
          critTermPoint[Re @ xSharp[a, b], sharpStyle]
          (* Re chops off small imaginary part in floating point arithmetic *)
        },
        (* One terminal point x_\[Natural] *)
        b == bNat[a],
        {
          critTermPoint[xNat[a], natStyle]
        },
        (* Zero critical terminal points *)
        True,
        {}
      ]
    ]
  , {{a, aInit, aIt}, aMin, aMax, Appearance -> "Open"}
  , {{b, bInit, bIt}, bMin, bMax, Appearance -> "Open"}
  , Button[resetLabel[aIt], a = aInit]
  , Button[resetLabel[bIt], b = bInit]
  ]
]


(* ::Subsection:: *)
(*Constructing an asymmetric domain*)


Module[
 {a, b,
  includeYReflection,
  xFl, xSh, sMax,
  xMin, xMax, yMax,
  mar, xMinMar, xMaxMar, yMaxMar,
  xInitGenericList, xyGenericList,
  xyInitInflectionList, xyInflectionList,
  xInitConvexList, sConvexLower, xyConvexList,
  emptyFrame,
  nonViableDomain,
  straightContour,
  lineOfSymmetry,
  inflectionFrontiers
 },
  (*
    ------------------------------------------------
    Definitions
    ------------------------------------------------
   *)
  (* A and B *)
  a = aAsymm;
  b = bAsymm;
  (* Include reflection in y (across x-axis) *)
  includeYReflection = {#, # * {1, -1}} &;
  (* Useful constants *)
  xFl = xFlat[a, b];
  xSh = xSharp[a, b];
  sMax = 6;
  (* Plot range *)
  xMin = Way[0, xFl] // N;
  xMax = 2;
  yMax = 2/3 (xMax - xMin);
  (* Plot range with margin *)
  mar = 1/5 xMax;
  xMinMar = xMin - mar;
  xMaxMar = xMax + mar;
  yMaxMar = yMax + mar;
  (* Generic traced boundaries from axis *)
  xInitGenericList = Subdivide[xFl, xSh, 4];
  xyGenericList = Table[
    With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
      NDSolveValue[
        {
          xyTraSystem[a, b],
          x[0] == xInit, y[0] == 0,
          WhenEvent[
            {
              x[s] < xMinMar,
              x[s] > xMaxMar,
              Abs @ y[s] > xMaxMar,
              vi[a, b][x[s], y[s]] < 0,
              tKnown[b][x[s], y[s]] < 0
            },
            "StopIntegration"
          ]
        }, {x, y}, {s, -sMax, sMax},
        NoExtrapolation
      ]
    ]
  , {xInit, xInitGenericList}];
  (* Inflection frontiers for the lower branch *)
  xyInitInflectionList = {
    {xInflAxis, 0},
    {xStraight, yInflStraight}
  };
  xyInflectionList = Table[
    With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
      NDSolveValue[
        {
          curContourSystem[a, b],
          x[0] == xyInit[[1]], y[0] == xyInit[[2]],
          WhenEvent[
            {
              x[s] < xMinMar,
              x[s] > xMaxMar,
              Abs @ y[s] > xMaxMar,
              vi[a, b][x[s], y[s]] < 0
            },
            "StopIntegration"
          ]
        }, {x, y}, {s, -sMax, sMax},
        NoExtrapolation
      ]
    ]
  , {xyInit, xyInitInflectionList}];
  (* Convex portions of traced boundaries *)
  xInitConvexList = Subdivide[xInflStraight, xInflAxis, 4];
  xyConvexList = Table[
    With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
      (*
        NOTE: WhenEvent detection doesn't work at the initial point,
        hence the If hack.
       *)
      sConvexLower = If[xInit == Last[xInitConvexList], 0, -sMax];
      NDSolveValue[
        {
          xyTraSystem[a, b],
          x[0] == xInit, y[0] == 0,
          WhenEvent[
            curTra[a, b, True][x[s], y[s]] < 0,
            "StopIntegration"
          ]
        }, {x, y}, {s, sConvexLower, sMax},
        NoExtrapolation
      ]
    ]
  , {xInit, xInitConvexList}];
  (*
    ------------------------------------------------
    Re-used plots
    ------------------------------------------------
   *)
  (* Empty frame *)
  emptyFrame =
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 360,
      PlotLabel -> BoxedLabel[
        Row[
          {aIt == a, bIt == N[b]},
          ","
        ]
      ]
    ];
  (* Non-viable domain *)
  nonViableDomain =
    RegionPlot[vi[a, b][x, y] < 0,
      {x, xMinMar, xMaxMar}, {y, -yMaxMar, yMaxMar},
      BoundaryStyle -> termStyle,
      PlotPoints -> 70,
      PlotStyle -> nonStyle
    ];
  (* Straight contour *)
  straightContour =
    Graphics @ {straightStyle,
      Line @ {{xStraight, -yMaxMar}, {xStraight, yMaxMar}}
    };
  (* Line of symmetry (y == 0) *)
  lineOfSymmetry =
    Graphics @ {Directive[Thin, Orange],
      Line @ N @ {{xFlat[a, b], 0}, {xSharp[a, b], 0}}
    };
  (* Inflection frontiers for the lower branch *)
  inflectionFrontiers[
    opt : OptionsPattern @ {"Mirror" -> False}
  ] :=
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // If[OptionValue["Mirror"], includeYReflection, Identity]
          // Evaluate,
        {s, DomainStart[xy], DomainEnd[xy]},
        PlotStyle -> inflStyle
      ]
    , {xy, xyInflectionList}];
  (*
    ------------------------------------------------
    Plots
    ------------------------------------------------
   *)
  {
    (* Generic traced boundaries from axis *)
    Show[
      emptyFrame,
      nonViableDomain,
      straightContour,
      lineOfSymmetry,
      Table[
        ParametricPlot[
          xy[s]
            // Through
            // includeYReflection
            // Evaluate,
          {s, DomainStart[xy], DomainEnd[xy]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {xy, xyGenericList}],
      inflectionFrontiers[]
    ]
      // Ex["cosine_general-asymmetric-traced-generic.pdf"],
    (* Convex portions of traced boundaries *)
    Show[
      emptyFrame,
      nonViableDomain,
      straightContour,
      lineOfSymmetry,
      Table[
        ParametricPlot[
          xy[s]
            // Through
            // includeYReflection
            // Evaluate,
          {s, DomainStart[xy], DomainEnd[xy]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {xy, xyConvexList}],
      inflectionFrontiers["Mirror" -> True]
    ]
      // Ex["cosine_general-asymmetric-traced-convex_portions.pdf"],
    (* Convex asymmetric domain boundaries *)
    Show[
      emptyFrame,
      nonViableDomain,
      straightContour,
      lineOfSymmetry,
      ParametricPlot[
        Table[
          xyTraAsymm[id][s] // Through
        , {id, {"upper", "lower"}}] // Evaluate,
        {s, -sMax, sMax},
        PlotStyle -> convexStyle
      ],
      inflectionFrontiers["Mirror" -> True]
    ]
      // Ex["cosine_general-asymmetric-traced-convex_domain.pdf"]
  }
]


(* ::Section:: *)
(*Numerical verification (B = 1) plots*)


(* ::Subsection:: *)
(*Finite element mesh*)


Table[
  Module[
   {source,
    a, mesh, prRad, prBath,
    bCoords, bCoordsRad, bCoordsBath,
    dest
   },
    (* Import mesh *)
    source = FString[
      "cosine_simple-verification-mesh-{aNamesSimpConvex[a]}.txt"
    ];
    {a, mesh, prRad, prBath} = Import[source] // Uncompress;
    (* Boundary coordinates *)
    bCoords = Part[
      mesh["Coordinates"],
      List @@ First @ mesh["BoundaryElements"] // Flatten // DeleteDuplicates
    ];
    bCoordsRad = Select[bCoords, prRad @@ # &];
    bCoordsBath = Select[bCoords, prBath @@ # &];
    (* Export plot *)
    dest = StringReplace[source, ".txt" -> ".pdf"];
    Show[
      mesh["Wireframe"],
      ListPlot[{bCoordsRad, bCoordsBath},
        PlotStyle -> (
          Directive[#, pointStyle, Opacity[0.7]] &
            /@ {Blue, Red}
        )
      ],
      ImageSize -> 240
    ] // Ex[dest]
  ]
, {a, aValuesSimpConvex}]


(* ::Subsection:: *)
(*Numerical solution*)


Table[
  Module[{source, tSol, mesh, dest},
    (* Import solution *)
    source = FString[
      "cosine_simple-verification-solution-{aNamesSimpConvex[a]}.txt"
    ];
    tSol = Import[source] // Uncompress;
    mesh = tSol["ElementMesh"];
    (* Plot *)
    dest = source // StringReplace[".txt" -> ".png"];
    With[{x = \[FormalX], y = \[FormalY]},
      Plot3D[tSol[x, y], Element[{x, y}, mesh],
        AxesLabel -> Italicise /@ {"x", "y", "T"},
        PlotLabel -> Column[
          {"Numerical solution", aIt == N[a]},
          Center
        ],
        PlotRange -> Full,
        PlotOptions[Axes] // Evaluate
      ]
    ] // Ex[dest]
  ]
, {a, aValuesSimpConvex}]


(* ::Subsection:: *)
(*Relative error (3D)*)


Table[
  Module[{source, tSol, mesh, dest},
    (* Import solution *)
    source = FString[
      "cosine_simple-verification-solution-{aNamesSimpConvex[a]}.txt"
    ];
    tSol = Import[source] // Uncompress;
    mesh = tSol["ElementMesh"];
    (* Plot *)
    dest = FString[
      "cosine_simple-verification-rel_error-3d-{aNamesSimpConvex[a]}.png"
    ];
    With[{x = \[FormalX], y = \[FormalY]},
      Plot3D[
        tSol[x, y] / tKnown[1][x, y] - 1, Element[{x, y}, mesh],
        PlotLabel -> Column[
          {"Relative error of numerical solution", aIt == N[a]},
          Center
        ],
        PlotRange -> Full,
        PlotOptions[Axes] // Evaluate
      ]
    ] // Ex[dest]
  ]
, {a, aValuesSimpConvex}]


(* ::Subsection:: *)
(*Relative error (2D)*)


Table[
  Module[{source, tSol, mesh, dest},
    (* Import solution *)
    source = FString[
      "cosine_simple-verification-solution-{aNamesSimpConvex[a]}.txt"
    ];
    tSol = Import[source] // Uncompress;
    mesh = tSol["ElementMesh"];
    (* Plot *)
    (* Note the aspect ratio *)
    dest = FString[
      "cosine_simple-verification-rel_error-2d-{aNamesSimpConvex[a]}.png"
    ];
    With[{x = \[FormalX], y = \[FormalY]},
      Show[
        DensityPlot[
          tSol[x, y] / tKnown[1][x, y] - 1, Element[{x, y}, mesh],
          ColorFunction -> "Rainbow",
          ImageSize -> 320,
          PlotLabel -> Column[
            {"Relative error of numerical solution", aIt == N[a]},
            Center
          ],
          PlotRange -> Full,
          PlotLegends -> Automatic,
          PlotOptions[Axes] // Evaluate
        ],
        mesh @ "Wireframe"[
          "MeshElementStyle" -> EdgeForm[Opacity[0.1]]
        ]
      ]
    ] // Ex[dest]
  ]
, {a, aValuesSimpConvex}]


(* ::Section:: *)
(*Numerical verification (B arbitrary) plots*)


(* ::Subsection:: *)
(*Finite element mesh*)


Module[
 {source,
  mesh, prRad, prBath,
  bCoords, bCoordsRad, bCoordsBath,
  dest
 },
  (* Import mesh *)
  source = "cosine_general-verification-mesh-asymmetric.txt";
  {mesh, prRad, prBath} = Import[source] // Uncompress;
  (* Boundary coordinates *)
  bCoords = Part[
    mesh["Coordinates"],
    List @@ First @ mesh["BoundaryElements"] // Flatten // DeleteDuplicates
  ];
  bCoordsRad = Select[bCoords, prRad @@ # &];
  bCoordsBath = Select[bCoords, prBath @@ # &];
  (* Export plot *)
  dest = StringReplace[source, ".txt" -> ".pdf"];
  Show[
    mesh["Wireframe"],
    ListPlot[{bCoordsRad, bCoordsBath},
      PlotStyle -> (
        Directive[#, pointStyle, Opacity[0.7]] &
          /@ {Blue, Red}
      )
    ],
    ImageSize -> 240
  ] // Ex[dest]
]


(* ::Subsection:: *)
(*Numerical solution*)


Module[{a, b, source, tSol, mesh, dest},
  a = aAsymm;
  b = bAsymm;
  (* Import solution *)
  source = "cosine_general-verification-solution-asymmetric.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Plot *)
  dest = source // StringReplace[".txt" -> ".png"];
  With[{x = \[FormalX], y = \[FormalY]},
    Plot3D[tSol[x, y], Element[{x, y}, mesh],
      AxesLabel -> Italicise /@ {"x", "y", "T"},
      PlotLabel -> Column[
        {
          "Numerical solution",
          Row[{aIt == a, bIt == N[b]}, ","]
        },
        Center
      ],
      PlotRange -> Full,
      PlotOptions[Axes] // Evaluate
    ]
  ] // Ex[dest]
]


(* ::Subsection:: *)
(*Relative error (3D)*)


Module[{a, b, source, tSol, mesh, dest},
  a = aAsymm;
  b = bAsymm;
  (* Import solution *)
  source = "cosine_general-verification-solution-asymmetric.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Plot *)
  dest = "cosine_general-verification-rel_error-3d-asymmetric.png";
  With[{x = \[FormalX], y = \[FormalY]},
    Plot3D[tSol[x, y] / tKnown[b][x, y] - 1, Element[{x, y}, mesh],
      AxesLabel -> Italicise /@ {"x", "y", "T"},
      PlotLabel -> Column[
        {
          "Rel. error of numerical solution",
          Row[{aIt == a, bIt == N[b]}, ","]
        },
        Center
      ],
      PlotRange -> Full,
      PlotOptions[Axes] // Evaluate
    ]
  ] // Ex[dest]
]


(* ::Subsection:: *)
(*Relative error (2D)*)


Module[{a, b, source, tSol, mesh, dest},
  a = aAsymm;
  b = bAsymm;
  (* Import solution *)
  source = "cosine_general-verification-solution-asymmetric.txt";
  tSol = Import[source] // Uncompress;
  mesh = tSol["ElementMesh"];
  (* Plot *)
  dest = "cosine_general-verification-rel_error-2d-asymmetric.png";
  With[{x = \[FormalX], y = \[FormalY]},
    Show[
      DensityPlot[tSol[x, y] / tKnown[b][x, y] - 1, Element[{x, y}, mesh],
        ColorFunction -> "Rainbow",
        PlotLabel -> Column[
          {
            "Rel. error of numerical solution",
            Row[{aIt == a, bIt == N[b]}, ","]
          },
          Center
        ],
        PlotRange -> Full,
        PlotLegends -> Automatic,
        PlotOptions[Frame] // Evaluate
      ],
      mesh["Wireframe"],
      (* Prevent bunched-up x ticks *)
      PlotRange -> {{1., 1.8}, {Automatic, Automatic}}
    ]
  ] // Ex[dest]
]


(* ::Section:: *)
(*Figure: (Un)physical region (cosine-physical.pdf)*)


Module[
 {bStep, bValues, bMax,
  xMin, xMax, yMax, imageSize,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinCont, xMaxCont, yMaxCont,
  yMaxContStraight,
  numTo1, numBeyond1,
  plotList,
  textStyle, arrowStyle,
  parameterArrow,
  xGraphicsBTick
 },
  (* Values of B *)
  bStep = 3/10;
  bValues = {1 - bStep, 1, 1 + bStep};
  bMax = Max[bValues];
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 3/2;
  yMax = 2;
  imageSize = 240;
  (* Margin *)
  eps = 0.05;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin - eps;
  xMaxUnphys = SeekRoot[tKnown[bMax][#, yMax] &, {0, xStraight}] + eps;
  yMaxUnphys = yMax + eps;
  (* Plot range for contours *)
  xMinCont = xMin - eps;
  xMaxCont = xMax + eps;
  yMaxCont = yMax + eps;
  (* Plot range for straight contour *)
  yMaxContStraight = yMax + 2 eps;
  (* Number of contours *)
  numTo1 = 5;
  numBeyond1 = 2;
  (* List of plots *)
  plotList =
    Table[
      Show[
        EmptyFrame[{xMin, xMax}, {-yMax, yMax},
          Frame -> None,
          ImageSize -> imageSize,
          PlotRangePadding -> None
        ],
        (* Unphysical domain *)
        RegionPlot[
          tKnown[b][x, y] < 0,
          {x, xMin, xMax}, {y, -yMax, yMax},
          BoundaryStyle -> BoundaryTracingStyle["Unphysical"],
          PlotPoints -> 10,
          PlotStyle -> BoundaryTracingStyle["Unphysical"]
        ],
        (* Known solution contours *)
        ContourPlot[
          tKnown[b][x, y],
          {x, xMinCont, xMaxCont}, {y, -yMaxCont, yMaxCont},
          ContourLabels -> None,
          Contours -> numTo1 + numBeyond1,
          ContourShading -> None,
          ContourStyle -> BoundaryTracingStyle["BackgroundDarker"],
          PlotPoints -> 5,
          PlotRange -> {0, 1 + (1 + numBeyond1) / numTo1}
        ],
        (* Straight contour *)
        ParametricPlot[{xStraight, y},
          {y, -yMaxContStraight, yMaxContStraight},
          PlotRange -> Full,
          PlotStyle -> BoundaryTracingStyle["ContourSolid"]
        ]
      ]
    , {b, bValues}];
  (* Parameter (B) increase indicator arrow *)
  textStyle = Style[#, 16] & @* LaTeXStyle;
  arrowStyle = Directive[Thickness[0.005], Arrowheads[0.04]];
  parameterArrow =
    Show[
      (* B-axis *)
      Graphics @ {arrowStyle,
        Arrow @ {{0, 0}, {1, 0}}
      },
      Graphics @ {
        Text[
          Italicise["B"] // textStyle,
          {1, 0},
          {-2, 0}
        ]
      },
      (* B == 1 *)
      xGraphicsBTick = 0.5;
      Graphics @ {arrowStyle,
        Line @ {
          {xGraphicsBTick, 0},
          {xGraphicsBTick, -0.01}
        }
      },
      Graphics @ {
        Text[
          1 // textStyle
          , {xGraphicsBTick, 0}
          , {0, 1.3}
        ]
      },
      ImageSize -> 1.85 imageSize,
      PlotRange -> All
    ];
  (* Final figure *)
  Column[
    {
      GraphicsRow[
        plotList,
        Spacings -> {0.2 imageSize, 0}
      ],
      parameterArrow
    },
    Spacings -> 0
  ]
] // Ex["cosine-physical.pdf"]


(* ::Section:: *)
(*Figure: Simple case (un)physical region and (non-)viable domain (cosine_simple-physical-viable.pdf)*)


Module[
 {b,
  aStep, aValues,
  xMin, xMax, yMax, imageSize,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinCont, xMaxCont, yMaxCont,
  xMinViable, xMaxViable, yMaxViable,
  yMaxContStraight,
  numTo1, numBeyond1,
  plotList,
  textStyle, arrowStyle,
  parameterArrow,
  xGraphicsATick
 },
  (* Value of B *)
  b = 1;
  (* Values of A *)
  aStep = 5/10;
  aValues = {1 - aStep, 1, 1 + aStep};
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 3/2;
  yMax = 2;
  imageSize = 240;
  (* Margin *)
  eps = 0.05;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin - eps;
  xMaxUnphys = SeekRoot[tKnown[b][#, yMax] &, {0, xStraight}] + eps;
  yMaxUnphys = yMax + eps;
  (* Plot range for contours *)
  xMinCont = xMin - eps;
  xMaxCont = xMax + eps;
  yMaxCont = yMax + eps;
  (* Plot range for straight contour *)
  yMaxContStraight = yMax + eps;
  (* Number of contours *)
  numTo1 = 5;
  numBeyond1 = 2;
  (* List of plots *)
  plotList =
    Table[
      (* Plot range for viable domain *)
      xMinViable = x0Simp[a] - eps;
      xMaxViable = xMax + eps;
      yMaxViable = yMax + eps;
      (* Plot *)
      Show[
        EmptyFrame[{xMin, xMax}, {-yMax, yMax},
          Frame -> None,
          ImageSize -> imageSize,
          PlotRangePadding -> None
        ],
        (* Unphysical domain *)
        RegionPlot[
          tKnown[b][x, y] < 0,
          {x, xMinUnphys, xMaxUnphys}, {y, -yMaxUnphys, yMaxUnphys},
          BoundaryStyle -> BoundaryTracingStyle["Unphysical"],
          PlotPoints -> 12,
          PlotStyle -> BoundaryTracingStyle["Unphysical"]
        ],
        (* Known solution contours *)
        ContourPlot[
          tKnown[b][x, y],
          {x, xMinCont, xMaxCont}, {y, -yMaxCont, yMaxCont},
          ContourLabels -> None,
          Contours -> numTo1 + numBeyond1,
          ContourShading -> None,
          ContourStyle -> BoundaryTracingStyle["BackgroundDarker"],
          PlotPoints -> 5,
          PlotRange -> {0, 1 + (1 + numBeyond1) / numTo1}
        ],
        (* Straight contour *)
        ParametricPlot[{xStraight, y},
          {y, -yMaxContStraight, yMaxContStraight},
          PlotRange -> Full,
          PlotStyle -> BoundaryTracingStyle["ContourSolid"]
        ],
        (* Non-viable domain *)
        RegionPlot[vi[a, b][x, y] < 0 && tKnown[b][x, y] > 0,
          {x, xMinViable, xMaxViable}, {y, -yMaxViable, yMaxViable},
          BoundaryStyle -> BoundaryTracingStyle["Terminal"],
          PlotPoints -> 8,
          PlotStyle -> BoundaryTracingStyle["NonViable"]
        ]
      ]
    , {a, aValues}];
  (* Parameter (A) increase indicator arrow *)
  textStyle = Style[#, 16] & @* LaTeXStyle;
  arrowStyle = Directive[Thickness[0.005], Arrowheads[0.04]];
  parameterArrow =
    Show[
      (* A-axis *)
      Graphics @ {arrowStyle,
        Arrow @ {{0, 0}, {1, 0}}
      },
      Graphics @ {
        Text[
          Italicise["A"] // textStyle,
          {1, 0},
          {-2, 0}
        ]
      },
      (* A == 1 *)
      xGraphicsATick = 0.54;
      Graphics @ {arrowStyle,
        Line @ {
          {xGraphicsATick, 0},
          {xGraphicsATick, -0.01}
        }
      },
      Graphics @ {
        Text[
          1 // textStyle
          , {xGraphicsATick, 0}
          , {0, 1.3}
        ]
      },
      ImageSize -> 1.85 imageSize,
      PlotRange -> All
    ];
  (* Final figure *)
  Column[
    {
      GraphicsRow[
        plotList,
        Spacings -> {0.2 imageSize, 0}
      ],
      parameterArrow
    },
    Spacings -> 0
  ]
] // Ex["cosine_simple-physical-viable.pdf"]


(* ::Section:: *)
(*Figure: Simple case traced boundaries (cosine_simple-traced-boundaries.pdf)*)


Module[
 {a, b,
  xMin, xMax, yMax, imageSize,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinViable, xMaxViable, yMaxViable,
  yMaxContStraight
 },
  (* Values of A and B *)
  a = 1/2;
  b = 1;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 3/2;
  yMax = 2;
  imageSize = 240;
  (* Margin *)
  eps = 0.05;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin - eps;
  xMaxUnphys = SeekRoot[tKnown[b][#, yMax] &, {0, xStraight}] + eps;
  yMaxUnphys = yMax + eps;
  (* Plot range for viable domain *)
  xMinViable = x0Simp[a] - eps;
  xMaxViable = xMax + eps;
  yMaxViable = yMax + eps;
  (* Plot range for straight contour *)
  yMaxContStraight = yMax + eps;
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> imageSize,
      LabelStyle -> LatinModernLabelStyle[12],
      PlotRangePadding -> None
    ],
    (* Unphysical domain *)
    RegionPlot[
      tKnown[b][x, y] < 0,
      {x, xMinUnphys, xMaxUnphys}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> BoundaryTracingStyle["Unphysical"],
      PlotPoints -> 12,
      PlotStyle -> BoundaryTracingStyle["Unphysical"]
    ],
    (* Non-viable domain *)
    RegionPlot[vi[a, b][x, y] < 0 && tKnown[b][x, y] > 0,
      {x, xMinViable, xMaxViable}, {y, -yMaxViable, yMaxViable},
      BoundaryStyle -> None,
      PlotPoints -> 2,
      PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Traced boundaries *)
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // {#, {#[[1]], -#[[2]]}} &
          // Evaluate,
        {s, DomainStart[xy], DomainEnd[xy]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
    , {xy, xyTraSimpFigure}]
  ]
] // Ex["cosine_simple-traced-boundaries.pdf"]


(* ::Section:: *)
(*Figure: Simple case traced boundaries, patched spiky (cosine_simple-traced-boundaries-patched-spiky.pdf)*)


Module[
 {a, b,
  xMin, xMax, yMax, imageSize,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinViable, xMaxViable, yMaxViable,
  yMaxContStraight
 },
  (* Values of A and B *)
  a = 1/2;
  b = 1;
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 3/2;
  yMax = 2;
  imageSize = 240;
  (* Margin *)
  eps = 0.05;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin - eps;
  xMaxUnphys = SeekRoot[tKnown[b][#, yMax] &, {0, xStraight}] + eps;
  yMaxUnphys = yMax + eps;
  (* Plot range for viable domain *)
  xMinViable = x0Simp[a] - eps;
  xMaxViable = xMax + eps;
  yMaxViable = yMax + eps;
  (* Plot range for straight contour *)
  yMaxContStraight = yMax + eps;
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      Frame -> None,
      ImageSize -> imageSize,
      PlotRangePadding -> None
    ],
    (* Traced boundaries (background) *)
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // Evaluate,
        {s, DomainStart[xy], DomainEnd[xy]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["BackgroundDarker"]
      ]
    , {xy, Join[patchedBoundaryUpperList, patchedBoundaryLowerList]}],
    (* Unphysical domain *)
    RegionPlot[
      tKnown[b][x, y] < 0,
      {x, xMinUnphys, xMaxUnphys}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> BoundaryTracingStyle["Unphysical"],
      PlotPoints -> 12,
      PlotStyle -> BoundaryTracingStyle["Unphysical"]
    ],
    (* Non-viable domain *)
    RegionPlot[vi[a, b][x, y] < 0 && tKnown[b][x, y] > 0,
      {x, xMinViable, xMaxViable}, {y, -yMaxViable, yMaxViable},
      BoundaryStyle -> None,
      PlotPoints -> 2,
      PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Traced boundaries (lower) *)
    Table[
      ParametricPlot[
        patchedBoundaryLowerList[[i]][s]
          // Through
          // Evaluate,
        {s, patchedIntersectionLowerList[[i]], 0}
        , PlotPoints -> 3
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
    , {i, patchedCornerNum}],
    (* Traced boundaries (upper) *)
    Table[
      ParametricPlot[
        patchedBoundaryUpperList[[i]][s]
          // Through
          // Evaluate,
        {s, 0, patchedIntersectionUpperList[[i]]}
        , PlotPoints -> 3
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
    , {i, patchedCornerNum}],
    {}
  ]
] // Ex["cosine_simple-traced-boundaries-patched-spiky.pdf"]


(* ::Section:: *)
(*Figure: Simple case traced boundaries, patched smooth (cosine_simple-traced-boundaries-patched-smooth.pdf)*)


Module[
 {a, b,
  x0,
  xMin, xMax, yMax, imageSize,
  eps,
  xMinUnphys, xMaxUnphys, yMaxUnphys,
  xMinViable, xMaxViable, yMaxViable,
  yMaxContStraight,
  textStyle, textStyleBracket
 },
  (* Values of A and B *)
  a = 1/2;
  b = 1;
  (* Critical terminal x-coordinate x_0 *)
  x0 = x0Simp[a];
  (* Plot range *)
  xMin = 0;
  xMax = Pi/2 * 3/2;
  yMax = 2;
  imageSize = 240;
  (* Margin *)
  eps = 0.05;
  (* Plot range for unphysical domain *)
  xMinUnphys = xMin - eps;
  xMaxUnphys = SeekRoot[tKnown[b][#, yMax] &, {0, xStraight}] + eps;
  yMaxUnphys = yMax + eps;
  (* Plot range for viable domain *)
  xMinViable = x0Simp[a] - eps;
  xMaxViable = xMax + eps;
  yMaxViable = yMax + eps;
  (* Plot range for straight contour *)
  yMaxContStraight = yMax + eps;
  (* Text style *)
  textStyle = Style[#, 18] & @* LaTeXStyle;
  textStyleBracket = Style[#, Larger] &;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      Frame -> None,
      ImageSize -> imageSize,
      PlotRangePadding -> None
    ],
    (* Traced boundaries (background) *)
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // {#, {#[[1]], -#[[2]]}} &
          // Evaluate,
        {s, DomainStart[xy], 0}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["BackgroundDarker"]
      ]
    , {xy, {patchedBoundarySmooth}}],
    (* Unphysical domain *)
    RegionPlot[
      tKnown[b][x, y] < 0,
      {x, xMinUnphys, xMaxUnphys}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> BoundaryTracingStyle["Unphysical"],
      PlotPoints -> 12,
      PlotStyle -> BoundaryTracingStyle["Unphysical"]
    ],
    (* Non-viable domain *)
    RegionPlot[vi[a, b][x, y] < 0 && tKnown[b][x, y] > 0,
      {x, xMinViable, xMaxViable}, {y, -yMaxViable, yMaxViable},
      BoundaryStyle -> None,
      PlotPoints -> 2,
      PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Traced boundaries (smooth patching) *)
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // {#, {#[[1]], -#[[2]]}} &
          // Evaluate,
        {s, 0, DomainEnd[xy]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
    , {xy, {patchedBoundarySmooth}}],
    (* Critical terminal point (x_0, 0) *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ {x0, 0}
    },
    Graphics @ {
      Text[
        Row @ {
          "(" // textStyleBracket,
          "\[NegativeVeryThinSpace]",
          Subscript[Italicise["x"], 0],
          ",\[ThinSpace]",
          0,
          ")" // textStyleBracket
        },
        {x0, 0},
        {1.5, -0.25}
      ]
        // LaTeXStyle
        // Style[#, 18] &
    },
    {}
  ]
] // Ex["cosine_simple-traced-boundaries-patched-smooth.pdf"]


(* ::Section:: *)
(*Figure: Simple case terminal points (cosine_simple-terminal-points.pdf)*)


Module[
 {a, b,
  x0,
  eps,
  xMin, xMax, yMax, imageSize,
  xMinCont, xMaxCont, yMaxCont,
  xMinViable, xMaxViable, yMaxViable,
  contNum, tContMin, tContStep, tContValues,
  tContOrd, xOrdGuess, yOrdGuess, xOrd, yOrd,
  textStyle, textStyleBracket, textVerticalShift
 },
  (* Values of A and B *)
  a = 1/2;
  b = 1;
  (* Critical terminal x-coordinate x_0 *)
  x0 = x0Simp[a];
  (* Margin *)
  eps = 0.01;
  (* Plot range *)
  xMin = x0 - 0.5 (xStraight - x0);
  xMax = xStraight;
  yMax = 1;
  imageSize = 240;
  (* Plot range for contours *)
  xMinCont = x0 - eps;
  xMaxCont = xMax + eps;
  yMaxCont = yMax + eps;
  (* Contours *)
  contNum = 4;
  tContMin = tKnown[b][x0, 0];
  tContStep = 0.03;
  tContValues = tContMin + tContStep (Range[contNum] - 1);
  (* Plot range for viable domain *)
  xMinViable = x0 - eps;
  xMaxViable = xMax + eps;
  yMaxViable = yMax + eps;
  (* Determine ordinary terminal point *)
  tContOrd = tContMin + tContStep;
  xOrdGuess = Way[x0, xMax];
  yOrdGuess = Way[0, yMax];
  {xOrd, yOrd} = FindRoot[
    Function[{x, y},
      {tKnown[b][x, y] - tContOrd, vi[a, b][x, y]}
    ],
    {{xOrdGuess}, {yOrdGuess}}
  ];
  (* Text style *)
  textStyle = Style[#, 18] & @* LaTeXStyle;
  textStyleBracket = Style[#, Larger] &;
  textVerticalShift = -0.25;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      AspectRatio -> 2,
      Frame -> None,
      ImageSize -> imageSize,
      PlotRangePadding -> None
    ],
    (* Known solution contours *)
    ContourPlot[
      tKnown[b][x, y],
      {x, xMinCont, xMaxCont}, {y, -yMaxCont, yMaxCont},
      ContourLabels -> None,
      Contours -> tContValues,
      ContourShading -> None,
      ContourStyle -> BoundaryTracingStyle["BackgroundDarker"],
      PlotPoints -> 8
    ],
    (* Non-viable domain *)
    RegionPlot[vi[a, b][x, y] < 0 && tKnown[b][x, y] > 0,
      {x, xMinViable, xMaxViable}, {y, -yMaxViable, yMaxViable},
      BoundaryStyle -> BoundaryTracingStyle["Terminal"],
      PlotPoints -> 9,
      PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Ordinary terminal point (x_ord, y_ord) *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ {xOrd, yOrd}
    },
    Graphics @ {
      Text[
        "ordinary",
        {xOrd, yOrd},
        {-1.4, textVerticalShift}
      ] // textStyle,
      {}
    },
    (* Critical terminal point (x_0, 0) *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ {x0, 0}
    },
    Graphics @ {
      Text[
        Row @ {
          "(" // textStyleBracket,
          "\[NegativeVeryThinSpace]",
          Subscript[Italicise["x"], 0],
          ",\[ThinSpace]",
          0,
          ")" // textStyleBracket
        },
        {x0, 0},
        {1.5, textVerticalShift}
      ] // textStyle,
      Text[
        "critical",
        {x0, 0},
        {-1.45, textVerticalShift}
      ] // textStyle,
      {}
    },
    {}
  ]
] // Ex["cosine_simple-terminal-points.pdf"]


(* ::Section:: *)
(*Figure: Simple case candidate domains (cosine_simple-candidate-domains.pdf)*)


Module[
  {
    b,
    numA,
    aStep, aValues,
    aMin, xRadAMin, yEndAMin,
    eps,
    yMaxFrame, xMinFrame, xMaxFrame,
    imageSize,
    plotList,
    xRad, yEnd,
    xRadEvenExtension,
    yInfl, xInfl,
    textStyle, arrowStyle,
    parameterArrow,
    xGraphicsAInfl, xGraphicsA1,
    dummyForTrailingCommas
   },
  (* Value of B *)
  b = 1;
  (* Number of A from A_i to 1 inclusive *)
  numA = 4;
  (* Values of A *)
  aStep = (1 - aInflSimp) / (numA - 1);
  aValues = aInflSimp + aStep * Range[-2, numA - 1];
  (* Plot range *)
  aMin = Min[aValues];
  xRadAMin = xTraCandSimp[aMin, True];
  yEndAMin = DomainEnd[xRadAMin];
  eps = 0.02;
  yMaxFrame = yEndAMin;
  xMinFrame = (1 - eps) x0Simp[aMin];
  xMaxFrame = (1 + eps) xStraight;
  (* List of plots *)
  imageSize = 96;
  plotList =
    Table[
      (* Radiation boundary x == x(y) for convex domain *)
      xRad = xTraCandSimp[a, True];
      yEnd = DomainEnd[xRad];
      xRadEvenExtension = Function[{y}, xRad[Abs @ y] // Evaluate];
      (* Point of inflection *)
      Which[
        a < aInflSimp,
          yInfl = SeekRoot[
            curTra[a, b][xRad[#], #] &,
            {DomainStart[xRad], DomainEnd[xRad]}
          ],
        a == aInflSimp,
          yInfl = DomainEnd[xRad],
        True,
          yInfl = Indeterminate
      ];
      xInfl = xRad[yInfl];
      (* Plot *)
      Show[
        EmptyFrame[
          {xMinFrame, xMaxFrame}, {-yMaxFrame, yMaxFrame}
          , Frame -> None
          , ImageSize -> 4 imageSize
            (* NOTE: GraphicsRow makes the plot smaller *)
          , PlotRangePadding -> None
        ],
        (* Convex domain *)
        ParametricPlot[
          {
            {xRadEvenExtension[y], y},
            {xStraight, y}
          }
          , {y, -yEnd, yEnd}
          , PlotPoints -> 2
          , PlotStyle -> BoundaryTracingStyle /@ {"Traced", "Contour"}
        ],
        (* Point of inflection *)
        If[NumericQ[yInfl],
          Graphics @ {GeneralStyle["Point"],
            Point @ {{xInfl, yInfl}, {xInfl, -yInfl}}
          },
          {}
        ],
        {}
      ]
      , {a, aValues}
    ];
  (* Parameter (A) increase indicator arrow *)
  textStyle = Style[#, 20] & @* LaTeXStyle;
  arrowStyle = Directive[Thickness[0.005], Arrowheads[0.04]];
  parameterArrow =
    Show[
      (* A-axis *)
      Graphics @ {arrowStyle,
        Arrow @ {{0, 0}, {1, 0}}
      },
      Graphics @ {
        Text[
          Italicise["A"] // textStyle
          , {1, 0}
          , {-2, 0}
        ]
      },
      (* A == A_infl *)
      xGraphicsAInfl = 0.41;
      Graphics @ {arrowStyle,
        Line @ {
          {xGraphicsAInfl, 0},
          {xGraphicsAInfl, -0.01}
        }
      },
      Graphics @ {
        Text[
          Subscript[Italicise["A"], "i"] // textStyle
          , {xGraphicsAInfl, 0}
          , {0, 1.3}
        ]
      },
      (* A == 1 *)
      xGraphicsA1 = 0.863;
      Graphics @ {arrowStyle,
        Line @ {
          {xGraphicsA1, 0},
          {xGraphicsA1, -0.01}
        }
      },
      Graphics @ {
        Text[
          1 // textStyle
          , {xGraphicsA1, 0}
          , {0, 1.3}
        ]
      },
      (* Convexity-interval labels *)
      Graphics @ {
        Text[
          "non\[Hyphen]convex" // textStyle
          (* NOTE: use '\\[Hyphen]' to prevent parsing as minus sign *)
          , {Way[0, xGraphicsAInfl], 0}
          , {0, -1.3}
        ],
        Text[
          "convex" // textStyle
          , {Way[xGraphicsAInfl, xGraphicsA1], 0}
          , {0, -1.3}
        ],
        {}
      },
      {}
      , ImageSize -> 4.2 imageSize
      , PlotRange -> All
    ];
  (* Final figure *)
  Column[
    {
      Row @ {
        GraphicsRow[
          plotList
          , Spacings -> {3.2 imageSize, 0}
        ],
        (* Adjust horizontal position *)
        Row @ {
          Graphics[ImageSize -> 0.5 imageSize]
        }
      },
      parameterArrow
    }
    , Center
    , Spacings -> 0
  ]
] // Ex["cosine_simple-candidate-domains.pdf"]


(* ::Section:: *)
(*Figure: General case (un)physical region and (non-)viable domain (cosine_general-physical-viable.pdf)*)


Module[
  {
    a, bN, bValues,
    xMin, xMax, yMax, imageSize,
    eps,
    plotList,
    xMinUnphys, xMaxUnphys, yMaxUnphys,
    plotPointsUnphys,
    xMinViable, xMaxViable, yMaxViable,
    plotPointsViable,
    textStyle, textStyleBracket, arrowStyle,
    parameterArrow,
    xGraphicsBNat, xGraphicsB1,
    dummyForTrailingCommas
  },
  (* Value of A *)
  a = 1;
  (* Value of B at gentle-to-fair transition *)
  bN = bNat[a];
  (* Values of B *)
  bValues = {
    (* Case 1 *) 0.7 bN,
    (* Case 2 *) bN,
    (* Case 3 *) Way[bN, 1, 0.03],
    (* Case 4 *) 1,
    (* Case 5 *) 2,
    Nothing
  };
  (* Plot range *)
  xMin = 0;
  xMax = 2;
  yMax = 3;
  imageSize = 240;
  (* Margin *)
  eps = 0.1;
  (* List of plots *)
  plotList =
    Table[
      (* Plot range for unphysical domain *)
      xMinUnphys = xMin - eps;
      xMaxUnphys = SeekRoot[tKnown[b][#, yMax] &, {0, xStraight}] + eps;
      yMaxUnphys = yMax + eps;
      plotPointsUnphys = If[b == 1, 9, 5];
      (* Plot range for viable domain *)
      xMinViable = If[b < 1, xMin - eps, xSharp[a, b] - eps];
      xMaxViable = xMax + eps;
      yMaxViable = yMax + eps;
      plotPointsViable = If[b < 1, 7, 8];
      (* Plot *)
      Show[
        EmptyFrame[{xMin, xMax}, {-yMax, yMax}
          , Frame -> None
          , ImageSize -> imageSize
          , PlotRangePadding -> None
        ],
        (* Unphysical domain *)
        RegionPlot[
          tKnown[b][x, y] < 0
          , {x, xMinUnphys, xMaxUnphys}
          , {y, -yMaxUnphys, yMaxUnphys}
          , BoundaryStyle -> BoundaryTracingStyle["Unphysical"]
          , PlotPoints -> plotPointsUnphys
          , PlotStyle -> BoundaryTracingStyle["Unphysical"]
        ],
        (* Non-viable domain *)
        RegionPlot[
          vi[a, b][x, y] < 0 && tKnown[b][x, y] > 0
          , {x, xMinViable, xMaxViable}
          , {y, -yMaxViable, yMaxViable}
          , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
          , PlotPoints -> plotPointsViable
          , PlotStyle -> BoundaryTracingStyle["NonViable"]
        ],
        (* Cirtical terminal points *)
        Which[
          b == bN,
          {
            Graphics @ {GeneralStyle["Point"],
              Point @ {N @ xNat[a], 0}
            },
            Graphics @ {
              Text[
                Subscript[Italicise["x"], "\[Natural]"] // textStyle
                , {N @ xNat[a], 0}
                , {0, -1.5}
             ]
            },
            {}
          },
          bN < b < 1,
          {
            Graphics @ {GeneralStyle["Point"],
              Point @ {N @ xFlat[a, b], 0}
            },
            Graphics @ {
              Text[
                Subscript[Italicise["x"], "\[Flat]"] // textStyle
                , {N @ xFlat[a, b], 0}
                , {0, -1.5}
             ]
            },
            Graphics @ {GeneralStyle["Point"],
              Point @ {N @ xSharp[a, b], 0}
            },
            Graphics @ {
              Text[
                Subscript[Italicise["x"], "\[Sharp]"] // textStyle
                , {N @ xSharp[a, b], 0}
                , {0.7, -1.5}
             ]
            },
            {}
          },
          b == 1,
          {
            Graphics @ {LightGray, GeneralStyle["Translucent"], GeneralStyle["Point"],
              Point @ {N @ xFlat[a, b], 0}
            },
            Graphics @ {
              Text[
                Subscript[Italicise["x"], "\[Flat]"] // textStyle
                , {N @ xFlat[a, b], 0}
                , {-2.3, -0.3}
             ]
            },
            Graphics @ {GeneralStyle["Point"],
              Point @ {N @ xSharp[a, b], 0}
            },
            Graphics @ {
              Text[
                Subscript[Italicise["x"], "\[Sharp]"] // textStyle
                , {N @ xSharp[a, b], 0}
                , {1.5, -1.1}
             ]
            },
            {}
          },
          b > 1,
          {
            Graphics @ {White, GeneralStyle["Point"],
              Point @ {N @ xFlat[a, b], 0}
            },
            Graphics @ {White,
              Text[
                Subscript[Italicise["x"], "\[Flat]"] // textStyle
                , {N @ xFlat[a, b], 0}
                , {0, -1.2}
             ]
            },
            Graphics @ {GeneralStyle["Point"],
              Point @ {N @ xSharp[a, b], 0}
            },
            Graphics @ {
              Text[
                Subscript[Italicise["x"], "\[Sharp]"] // textStyle
                , {N @ xSharp[a, b], 0}
                , {1.5, -1.1}
             ]
            },
            {}
          },
          True, {}
        ],
        {}
      ]
      , {b, bValues}
    ];
  (* Parameter (B) increase indicator arrow *)
  textStyle = Style[#, 16] & @* LaTeXStyle;
  textStyleBracket = Style[#, 22] &;
  arrowStyle = Directive[Thickness[0.005], Arrowheads[0.04]];
  parameterArrow =
    Show[
      (* B-axis *)
      Graphics @ {arrowStyle,
        Arrow @ {{0, 0}, {1, 0}}
      },
      Graphics @ {
        Text[
          Italicise["B"] // textStyle,
          {1, 0},
          {-2, 0}
        ]
      },
      (* B == B_nat *)
      xGraphicsBNat = 0.288;
      Graphics @ {arrowStyle,
        Line @ {
          {xGraphicsBNat, 0},
          {xGraphicsBNat, -0.01}
        }
      },
      Graphics @ {
        Text[
          Row @ {
            Subscript[Italicise["B"], "\[Natural]"],
            "(" // textStyleBracket,
            "\[NegativeVeryThinSpace]",
            Italicise["A"],
            ")" // textStyleBracket,
            Nothing
          }
            // textStyle
          , {xGraphicsBNat, 0}
          , {0, 1.4}
        ]
      },
      (* B == 1 *)
      xGraphicsB1 = 0.725;
      Graphics @ {arrowStyle,
        Line @ {
          {xGraphicsB1, 0},
          {xGraphicsB1, -0.01}
        }
      },
      Graphics @ {
        Text[
          1 // textStyle
          , {xGraphicsB1, 0}
          , {0, 1.4}
        ]
      },
      (* Regime labels *)
      Graphics @ {
        Text[
          "gentle" // textStyle
          , {Way[0, xGraphicsBNat, 0.25], 0}
          , {0, -1.3}
        ],
        Text[
          "fair" // textStyle
          , {Way[xGraphicsBNat, xGraphicsB1, 0.48], 0}
          , {0, -1.3}
        ],
        Text[
          "steep" // textStyle
          , {Way[xGraphicsB1, 1, 0.69], 0}
          , {0, -1.3}
        ],
        {}
      },
      {}
      , ImageSize -> 2.07 imageSize
      , PlotRange -> All
    ];
  (* Final figure *)
  Column[
    {
      GraphicsRow[
        plotList,
        Spacings -> {0.3 imageSize, 0}
      ],
      parameterArrow
    },
    Spacings -> 0
  ]
] // Ex["cosine_general-physical-viable.pdf"]


(* ::Section:: *)
(*Figure: General case critical terminal points (cosine_general-critical.pdf)*)


Module[
  {
    a,
    bN, xN,
    bMin, bMax,
    xMin, xMax, xMaxFilling,
    xFillingFunction,
    dummyForTrailingCommas
  },
  (* Value of A *)
  a = 1;
  (* Critical values B_\[Natural] and x_\[Natural] *)
  (* (N needed otherwise Epilog doesn't work) *)
  bN = bNat[a] // N;
  xN = xNat[a] // N;
  (* Plot range *)
  bMin = 0;
  bMax = 1.7;
  xMin = 0;
  xMax = 2;
  xMaxFilling = Way[Pi/2, xMax, 2/3];
  (* Filling function for B < B_\[Natural] *)
  xFillingFunction[b_] :=
    Piecewise @ {
      {xMin, b < bN},
      {Indeterminate, True}
    };
  (* Plot *)
  Plot[
    {xSharp[a, b], xFlat[a, b], xFillingFunction[b]}
    , {b, bMin, bMax}
    , AspectRatio -> 1 / 1.5
    , AxesLabel -> {bIt, xIt}
    , Epilog -> {
        (* Guiding lines for critical values *)
        BoundaryTracingStyle["Contour"],
        Line @ {{0, xN}, {bN, xN}, {bN, 0}},
        (* Points for regime transitions *)
        GeneralStyle["Point"],
        Point @ {{bN, xN}, {1, 0}},
        {}
      }
    , Filling -> {1 -> xMaxFilling, 2 -> xMin, 3 -> xMaxFilling}
    , FillingStyle -> BoundaryTracingStyle["NonViable"]
    , ImageSize -> 240
    , LabelStyle -> LatinModernLabelStyle[16]
    , PlotLegends -> Placed[
        {
          xIt == Subscript[xIt, "\[Sharp]"],
          xIt == Subscript[xIt, "\[Flat]"],
          Nothing
        },
        After
      ]
    , PlotRange -> {xMin, xMax}
    , PlotRangeClipping -> False
    , PlotStyle -> {Black, Directive[Black, Dashed], None}
    , Ticks -> {
        {0, {bN, Subscript[bIt, "\[Natural]"]}, 1},
        {0, {xN, Subscript[xIt, "\[Natural]"]}}
      }
  ]
] // Ex["cosine_general-critical.pdf"]



