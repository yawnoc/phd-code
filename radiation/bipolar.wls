#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*See r3 (manuscripts/radiation-3-bipolar.pdf).*)


(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< Curvilinear`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "bipolar"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Critical terminal points along u = \[Pi]*)


(* ::Subsubsection:: *)
(*(v_\[Natural]\[Pi], A_\[Natural]\[Pi])*)
(*(Warm-to-hot transition hyperbolic coordinate and dimensionless group)*)


(* ::Text:: *)
(*See Section 4.2.1 of r3 (Page r3-8), particularly (r3.43) and (r3.44).*)


vNatPi =
  InverseFunction[
    Function[{v},
      ConditionalExpression[
        Cosh[v] + 1 - v Sinh[v] / 4,
        v > 0
      ]
    ]
  ][0];


aNatPi = Function[{v}, 4 v^3 / Sinh[v]] @ vNatPi;


{vNatPi, aNatPi} // N


(* ::Subsubsection:: *)
(*v_\[Flat]\[Pi] and v_\[Sharp]\[Pi]*)


(* ::Text:: *)
(*See Sections 4.1.1 and 4.1.2 of r3 (Page r3-6).*)


(* The argument a to vFlatPi and vSharpPi must be EXACT. *)


vFlatPi[a_] /; a <= aNatPi :=
  InverseFunction[
    Function[{v},
      ConditionalExpression[
        Cosh[v] + 1 - v^4 / a,
        0 < v <= vNatPi
      ] // Evaluate
    ]
  ][0] // Evaluate;


vSharpPi[a_] /; a <= aNatPi :=
  InverseFunction[
    Function[{v},
      ConditionalExpression[
        Cosh[v] + 1 - v^4 / a,
        v >= vNatPi
      ] // Evaluate
    ]
  ][0] // Evaluate;


(* Check warm-to-hot transition (u = \[Pi]) *)
vFlatPi[aNatPi] == vNatPi == vSharpPi[aNatPi]


(* ::Subsection:: *)
(*Critical terminal points along u = 0*)


(* ::Subsubsection:: *)
(*(v_\[Natural]0, A_\[Natural]0)*)
(*(Cold-to-warm transition hyperbolic coordinate and dimensionless group)*)


(* ::Text:: *)
(*See Section 4.2.2 of r3 (Page r3-8), particularly (r3.49) and (r3.50).*)


vNat0 =
  InverseFunction[
    Function[{v},
      ConditionalExpression[
        Cosh[v] - 1 - v Sinh[v] / 4,
        v > 0
      ]
    ]
  ][0];


aNat0 = Function[{v}, 4 v^3 / Sinh[v]] @ vNat0;


{vNat0, aNat0} // N


(* ::Subsubsection:: *)
(*v_\[Flat]0 and v_\[Sharp]0*)


(* ::Text:: *)
(*See Sections 4.1.1 to 4.1.4 of r3 (Pages r3-6 and r3-7).*)


(* The argument a to vFlat0 and vSharp0 must be EXACT. *)


vFlat0[a_] /; a <= aNat0 :=
  InverseFunction[
    Function[{v},
      ConditionalExpression[
        Cosh[v] - 1 - v^4 / a,
        0 < v <= vNat0
      ] // Evaluate
    ]
  ][0] // Evaluate;


vSharp0[a_] /; a <= aNat0 :=
  InverseFunction[
    Function[{v},
      ConditionalExpression[
        Cosh[v] - 1 - v^4 / a,
        v >= vNat0
      ] // Evaluate
    ]
  ][0] // Evaluate;


(* Check cold-to-warm transition (u = 0) *)
vFlat0[aNat0] == vNat0 == vSharp0[aNat0]


(* ::Subsection:: *)
(*\[CapitalPhi] (viability)*)


(* ::Text:: *)
(*See (r3.39) (Page r3-5).*)


vi[a_] := Function[{u, v},
  (Cosh[v] - Cos[u])^2 - v^8 / a^2
];


(* ::Subsection:: *)
(*\[Omega] = \[Omega](v_\[Sharp]0)*)


(* ::Text:: *)
(*See (r3.86) (Page r3-14).*)


omega[v_] := (Cosh[v] - 1) / (v Sinh[v]);


(* ::Subsection:: *)
(*Representative values of A*)


aHot = 9;
aWarmHot = aNatPi;
aWarm = 9 + 1/2;
aColdWarm = aNat0;
aCold = 10;
aHot < aWarmHot < aWarm < aColdWarm < aCold


(* ::Subsection:: *)
(*Traced boundary derivative dv/du*)


(* ::Text:: *)
(*See (r3.53) (Page r3-9).*)


vTraDer[a_] := Function[{u, v},
  (a / v^4) Sqrt @ vi[a][u, v] // Evaluate
] // Evaluate;


(* ::Subsection:: *)
(*Traced boundary curvature*)


(* ::Text:: *)
(*See (r3.63) (Page r3-10).*)


curTra[a_] := Function[{u, v},
  Module[{d, dd, h, cur, vDer, vDer2},
    (* Abbreviation for u-derivative *)
    d = Dt[#, u, Constants -> a] &;
    (* Abbreviations *)
    h = HBipolar[u, v];
    dd = Sinh[v] - d[v] Sin[u];
    (* Curvature *)
    cur = h^3 (dd (1 + d[v]^2) + d @ d[v] / h);
    (* v' and v'' *)
    vDer = vTraDer[a][u, v];
    vDer2 = d[vDer];
    (* Curvature evaluated *)
    cur
      /. {d @ d[v] -> vDer2}
      /. {d[v] -> vDer}
  ] // Evaluate
] // Evaluate;


(* ::Subsubsection:: *)
(*v_i (inflection hyperbolic coordinate along u = -\[Pi])*)


(* ::Text:: *)
(*See (r3.65) (Page r3-11).*)


vInfl = InverseFunction[
    Function[{v},
      ConditionalExpression[
        v/2 Tanh[v/2] - 1,
        v > 0
      ] // Evaluate
    ]
  ][0];


(* ::Subsubsection:: *)
(*Candidate boundary (at u = -\[Pi])*)


(* ::Text:: *)
(*See the bottom of Page r3-9.*)


(* The hyperbolic coordinate (v) at u == -\[Pi] *)
vTraCandCorn[a_] /; 0 < a <= aNat0 :=
  With[{v = \[FormalV]},
    NDSolveValue[
      {
        v'[u] == Re @ vTraDer[a][u, v[u]],
        v[0] == vSharp0[a],
        WhenEvent[vi[a][u, v[u]] < 0, "StopIntegration"]
      }, v[-Pi], {u, -Pi, 0},
      NoExtrapolation,
      PreciseOptions[]
    ]
  ];


(* ::Subsubsection:: *)
(*A_i (inflection dimensionless group)*)


(* ::Text:: *)
(*See (r3.68) and (r3.69) (Page r3-11).*)


(* Compute A_i using the bisection algorithm *)
aInfl = Module[{dest, a, aMin, aMax, num},
  (* (This is really slow (~30 sec), so compute once and store.) *)
  (* (Delete the file manually to compute from scratch.) *)
  dest = "bipolar-a-inflection.txt";
  If[FileExistsQ[dest],
    (* If already stored, import *)
    {a, num} = Import[dest] // Uncompress,
    (* Otherwise compute and export for next time *)
    (* (argument to vTraCandCorn must be EXACT) *)
    aMin = aNat0 - 4/1000;
    aMax = aNat0;
    {a, num} = SeekRootBisection[
      vTraCandCorn[#] - vInfl &,
      {aMin, aMax},
      10^-7,
      N -> False,
      "ReturnIterations" -> True
    ];
    {a, num} // Compress // Ex[dest]
  ];
  (* Print iterations used *)
  Print["Bisection algorithm: " <> ToString[num] <> " iterations"];
  (* Return value of A_i *)
  a
];
aInfl // N


(* ::Subsubsection:: *)
(*Asymmetry*)


(* ::Text:: *)
(*See (r3.92) (Page r3-16).*)


asymm[vSharp0_, vCorn_] :=
  Divide[
    1 - XBipolar[-Pi, vCorn],
    XBipolar[0, vSharp0] - 1
  ] // FullSimplify // Evaluate;


(* ::Subsection:: *)
(*Starting points for boundary tracing*)


(* ::Subsubsection:: *)
(*System of ODES for terminal curve \[CapitalPhi] = 0*)


(* ::Text:: *)
(*See (57) and (58) (Page 13 of manuscripts/boundary-tracing.pdf)*)
(*for arc length parametrisation of a contour.*)


viContourSystem[a_, sSign_: 1] :=
  With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
    Module[{h, fun, p, q, slope},
      (* Scale factor for bipolar coordinates *)
      h = HBipolar[u, v];
      (* Scalar function \[CapitalPhi] whose contours are sought *)
      fun = vi[a][u, v];
      (* Components of the gradient vector *)
      p = (1 / h) D[fun, u];
      q = (1 / h) D[fun, v];
      (* Magnitude of the gradient vector *)
      slope = Sqrt[p^2 + q^2];
      (* Return system of ODEs *)
      {
        u' == q / (h slope),
        v' == -p / (h slope)
      } /. {
        u' -> Sign[sSign] u'[s],
        v' -> Sign[sSign] v'[s],
        u -> u[s],
        v -> v[s],
        List -> Sequence
      }
    ]
  ] // Evaluate;


(* ::Subsubsection:: *)
(*Hot regime A < A_\[Natural]\[Pi]*)


(* Starting points which are hyperbolic critical terminal points *)
(* (the two critical terminal points along u == 0) *)
startUVHot["hyperbolic"] =
  Module[{a},
    a = aHot;
    {
      {0, vSharp0[a]},
      {0, vFlat0[a]}
    } // N // Rationalize[#, 0] &
  ];


(* Starting points along outer terminal curve *)
startUVHot["outer"] =
  With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
    Module[{a, nmax, smax, uvTerm},
      a = aHot;
      nmax = 16;
      (* Upper bound for arc length traversed *)
      smax = 2 Pi (XBipolar[0, vFlat0[a]] - 1);
      (* Outer terminal curve *)
      uvTerm =
        NDSolveValue[
          {
            viContourSystem[a, -1],
            u[0] == 0, v[0] == vFlat0[a],
            WhenEvent[u[s] == 2 Pi, "StopIntegration"]
          }, {u, v}, {s, 0, smax},
          NoExtrapolation
        ];
      (* Actual arc length traversed *)
      smax = DomainEnd[uvTerm];
      (* Starting points on the outer terminal curve *)
      Table[
        uvTerm[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[0, smax, nmax]}] // Most
    ]
  ];


(* Starting points along inner terminal curve *)
startUVHot["inner"] =
  With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
    Module[{a, nmax, smax, uvTerm},
      a = aHot;
      nmax = 16;
      (* Upper bound for arc length traversed *)
      smax = 2 Pi (XBipolar[0, vSharpPi[a]] - 1);
      (* Inner terminal curve *)
      uvTerm =
        NDSolveValue[
          {
            viContourSystem[a],
            u[0] == 0, v[0] == vSharp0[a],
            WhenEvent[u[s] == 2 Pi, "StopIntegration"]
          }, {u, v}, {s, 0, smax},
          NoExtrapolation
        ];
      (* Actual arc length traversed *)
      smax = DomainEnd[uvTerm];
      (* Starting points on the outer terminal curve *)
      Table[
        uvTerm[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[0, smax, nmax]}] // Most
    ]
  ];


(* Starting points along axis *)
startUVHot["axis"] =
  Module[{a, xOuter, xValues},
    a = aHot;
    xOuter = XBipolar[0, vFlat0[a]];
    xValues =
      Join[
        Subdivide[1, 17/16, 4] // Rest,
        {9/8, 5/4, 3/2, 2}
      ] xOuter;
    Table[{0, VBipolar[x, 0]}, {x, xValues}]
  ];


(* ::Subsubsection:: *)
(*Warm-to-hot transition A = A_\[Natural]\[Pi]*)


(* Starting points which are hyperbolic critical terminal points *)
(* (the two critical terminal points along u == 0) *)
startUVWarmHot["hyperbolic"] =
  Module[{a},
    a = aWarmHot;
    {
      {0, vSharp0[a]},
      {0, vFlat0[a]}
    } // N // Rationalize[#, 0] &
  ];


(* Starting points along outer terminal curve *)
startUVWarmHot["outer"] =
  With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
    Module[{a, nmax, smax, vPert, uvTerm},
      a = aWarmHot;
      nmax = 16;
      (* Upper bound for arc length traversed *)
      (* (between -smax and smax to avoid pincer at vNatPi) *)
      smax = Pi (XBipolar[0, vFlat0[a]] - 1);
      (* Prevent sharp turn at vNatPi *)
      vPert = 10^-6;
      (* Outer terminal curve *)
      uvTerm =
        NDSolveValue[
          {
            viContourSystem[a, -1],
            u[0] == 0, v[0] == vFlat0[a] - vPert,
            WhenEvent[Abs @ u[s] == Pi, "StopIntegration"]
          }, {u, v}, {s, -smax, smax},
          NoExtrapolation
        ];
      (* Actual arc length traversed *)
      smax = DomainEnd[uvTerm];
      (* Starting points on the outer terminal curve *)
      Table[
        uvTerm[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[-smax, smax, nmax]}]
        // Most
        // RotateLeft[#, nmax / 2] &
    ]
  ];


(* Starting points along inner terminal curve *)
startUVWarmHot["inner"] =
  With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
    Module[{a, nmax, smax, vPert, uvTerm},
      a = aWarmHot;
      nmax = 16;
      (* Upper bound for arc length traversed *)
      (* (between -smax and smax to avoid pincer at vNatPi) *)
      smax = Pi (XBipolar[0, vSharpPi[a]] - 1);
      (* Prevent sharp turn at vNatPi *)
      vPert = 10^-6;
      (* Inner terminal curve *)
      uvTerm =
        NDSolveValue[
          {
            viContourSystem[a],
            u[0] == 0, v[0] == vSharp0[a] + vPert,
            WhenEvent[Abs @ u[s] == Pi, "StopIntegration"]
          }, {u, v}, {s, -smax, smax},
          NoExtrapolation
        ];
      (* Actual arc length traversed *)
      smax = DomainEnd[uvTerm];
      (* Starting points on the outer terminal curve *)
      Table[
        uvTerm[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[-smax, smax, nmax]}]
        // Most
        // RotateLeft[#, nmax / 2] &
    ]
  ];


(* Starting points along axis *)
startUVWarmHot["axis"] =
  Module[{a, xOuter, xValues},
    a = aWarmHot;
    xOuter = XBipolar[0, vFlat0[a]];
    xValues =
      Join[
        Subdivide[1, 17/16, 4] // Rest,
        {9/8, 5/4, 3/2, 2}
      ] xOuter;
    Table[{0, VBipolar[x, 0]}, {x, xValues}]
  ];


(* ::Subsubsection:: *)
(*Warm regime A_\[Natural]\[Pi] < A < A_\[Natural]0*)


(* Starting points which are hyperbolic critical terminal points *)
(* (the two critical terminal points along u == 0) *)
startUVWarm["hyperbolic"] =
  Module[{a},
    a = aWarm;
    {
      {0, vSharp0[a]},
      {0, vFlat0[a]}
    } // N // Rationalize[#, 0] &
  ];


(* Starting points along terminal curve *)
startUVWarm["terminal"] =
  With[{u = \[FormalU], v = \[FormalV], s = \[FormalS]},
    Module[{a, nmax, smax, uvTerm},
      a = aWarm;
      nmax = 16;
      (* (Probable) upper bound for arc length traversed *)
      smax = 2 Pi (XBipolar[0, vFlat0[a]] - 1);
      (* Outer terminal curve *)
      uvTerm =
        NDSolveValue[
          {
            viContourSystem[a, -1],
            u[0] == 0, v[0] == vFlat0[a],
            WhenEvent[u[s] == 0 && u'[s] > 0, "StopIntegration"]
          }, {u, v}, {s, 0, smax},
          NoExtrapolation
        ];
      (* Actual arc length traversed *)
      smax = DomainEnd[uvTerm];
      (* Starting points on the outer terminal curve *)
      Table[
        uvTerm[s] // Through // Rationalize[#, 0] &
      , {s, Subdivide[0, smax, nmax]}] // Most
    ]
  ];


(* Starting points along axis *)
startUVWarm["axis"] =
  Module[{a, xOuter, xValues},
    a = aWarm;
    xOuter = XBipolar[0, vFlat0[a]];
    xValues =
      Join[
        Subdivide[1, 17/16, 4] // Rest,
        {9/8, 5/4, 3/2, 2}
      ] xOuter;
    Table[{0, VBipolar[x, 0]}, {x, xValues}]
  ];


(* ::Subsubsection:: *)
(*Cold-to-warm transition A = A_\[Natural]0*)


(* Starting point which is a hyperbolic critical terminal point *)
(* (the lone critical terminal point along u == 0) *)
startUVColdWarm["hyperbolic"] = {{0, vNat0}} // N // Rationalize[#, 0] &;


(* Starting points along terminal curve *)
(* (the terminal curve is degenerate: it is simply the point {0, vNat0}) *)


(* Starting points along axis *)
startUVColdWarm["axis"] =
  Module[{nmax, xTerminal, xValues},
    nmax = 16;
    xTerminal = XBipolar[0, vNat0];
    xValues = Subdivide[1, 2 xTerminal - 1, nmax];
    xValues = xValues // Rest;
    Table[{0, VBipolar[x, 0]}, {x, xValues}]
  ];


(* ::Subsubsection:: *)
(*Cold regime A > A_\[Natural]0*)


(* There are no critical terminal points; there is no terminal curve *)


(* Starting points along axis *)


startUVCold["axis"] =
  Module[{nmax, xValues},
    nmax = 16;
    xValues = Subdivide[1, 2, nmax] // Rest;
    Table[{0, VBipolar[x, 0]}, {x, xValues}]
  ];


(* ::Subsection:: *)
(*Traced boundary v = v(u)*)


(* ::Subsubsection:: *)
(*Hot regime A < A_\[Natural]\[Pi]*)


With[{v = \[FormalV]},
  Module[{a, uvInitList, viTol, uInit, vInit},
    a = aHot;
    Table[
      uvInitList = startUVHot[id];
      viTol = vi[a] @@@ uvInitList // Min;
      vTraHot[id] =
        Table[
          {uInit, vInit} = uvInit;
          NDSolveValue[
            {
              v'[u] == Re @ vTraDer[a][u, v[u]],
              v[uInit] == vInit,
              WhenEvent[vi[a][u, v[u]] < viTol, "StopIntegration"]
            }, v, {u, uInit - 2 Pi, uInit + 2 Pi},
            NoExtrapolation
          ]
        , {uvInit, uvInitList}];
    , {id, {"hyperbolic", "outer", "inner", "axis"}}]
  ];
];


(* ::Subsubsection:: *)
(*Warm-to-hot transition A = A_\[Natural]\[Pi]*)


With[{v = \[FormalV]},
  Module[{a, uvInitList, viTol, uInit, vInit},
    a = aWarmHot;
    Table[
      uvInitList = startUVWarmHot[id];
      viTol = vi[a] @@@ uvInitList // Min;
      vTraWarmHot[id] =
        Table[
          {uInit, vInit} = uvInit;
          NDSolveValue[
            {
              v'[u] == Re @ vTraDer[a][u, v[u]],
              v[uInit] == vInit,
              WhenEvent[vi[a][u, v[u]] < viTol, "StopIntegration"]
            }, v, {u, uInit - 2 Pi, uInit + 2 Pi},
            NoExtrapolation
          ]
        , {uvInit, uvInitList}];
    , {id, {"hyperbolic", "outer", "inner", "axis"}}]
  ];
];


(* ::Subsubsection:: *)
(*Warm regime A_\[Natural]\[Pi] < A < A_\[Natural]0*)


With[{v = \[FormalV]},
  Module[{a, uvInitList, viTol, uInit, vInit},
    a = aWarm;
    Table[
      uvInitList = startUVWarm[id];
      viTol = vi[a] @@@ uvInitList // Min;
      vTraWarm[id] =
        Table[
          {uInit, vInit} = uvInit;
          NDSolveValue[
            {
              v'[u] == Re @ vTraDer[a][u, v[u]],
              v[uInit] == vInit,
              WhenEvent[vi[a][u, v[u]] < viTol, "StopIntegration"]
            }, v, {u, uInit - 2 Pi, uInit + 2 Pi},
            NoExtrapolation
          ]
        , {uvInit, uvInitList}];
    , {id, {"hyperbolic", "terminal", "axis"}}]
  ];
];


(* ::Subsubsection:: *)
(*Cold-to-warm transition A = A_\[Natural]0*)


With[{v = \[FormalV]},
  Module[{a, uvInitList, viTol, uInit, vInit},
    a = aColdWarm;
    Table[
      uvInitList = startUVColdWarm[id];
      viTol = vi[a] @@@ uvInitList // Min;
      vTraColdWarm[id] =
        Table[
          {uInit, vInit} = uvInit;
          NDSolveValue[
            {
              v'[u] == Re @ vTraDer[a][u, v[u]],
              v[uInit] == vInit,
              WhenEvent[vi[a][u, v[u]] < viTol, "StopIntegration"]
            }, v, {u, uInit - 2 Pi, uInit + 2 Pi},
            NoExtrapolation
          ]
        , {uvInit, uvInitList}];
    , {id, {"hyperbolic", "axis"}}]
  ];
];


(* ::Subsubsection:: *)
(*Cold regime A > A_\[Natural]0*)


With[{v = \[FormalV]},
  Module[{a, uvInitList, viTol, uInit, vInit},
    a = aCold;
    Table[
      uvInitList = startUVCold[id];
      viTol = vi[a] @@@ uvInitList // Min;
      vTraCold[id] =
        Table[
          {uInit, vInit} = uvInit;
          NDSolveValue[
            {
              v'[u] == Re @ vTraDer[a][u, v[u]],
              v[uInit] == vInit,
              WhenEvent[vi[a][u, v[u]] < viTol, "StopIntegration"]
            }, v, {u, uInit - 2 Pi, uInit + 2 Pi},
            NoExtrapolation
          ]
        , {uvInit, uvInitList}];
    , {id, {"axis"}}]
  ];
];


(* ::Subsection:: *)
(*Italicised symbols*)


aIt = Italicised["A"];
uIt = Italicised["u"];
vIt = Italicised["v"];


(* ::Subsection:: *)
(*Options for exported GIFS*)


gifOpts = Sequence[
  AnimationRepetitions -> Infinity,
  "DisplayDurations" -> 0.5
];


(* ::Subsection:: *)
(*Bipolar coordinates*)


(* ::Subsubsection:: *)
(*Effectively v = Infinity*)


vInfinity = 10;


(* ::Subsubsection:: *)
(*Contours*)


defContStyle = LightGray;


uContour[u_, style_: defContStyle] :=
  ParametricPlot[
    XYBipolar[u, v], {v, -vInfinity, vInfinity},
    Axes -> None,
    Exclusions -> If[Mod[u, 2 Pi] == 0, 0, Automatic],
    PlotStyle -> style
  ];


vContour[v_, style_: defContStyle] :=
  ParametricPlot[
    XYBipolar[u, v], {u, 0, 2 Pi},
    Axes -> None,
    PlotStyle -> style
  ];


uContours = uContour /@ Most[Range[0, 360, 30] Degree];
vContours = vContour /@ {-2, -1, -1/2, 0, 1/2, 1, 2};


(* ::Subsubsection:: *)
(*Local orthonormal basis vectors*)


uVector[u_, v_, style_] :=
  Graphics @ {Directive[Thick, style],
    Arrow @ {
      XYBipolar[u, v],
      XYBipolar[u, v] + AUBipolar[u, v]
    }
  };


vVector[u_, v_, style_] :=
  Graphics @ {Directive[Thick, style],
    Arrow @ {
      XYBipolar[u, v],
      XYBipolar[u, v] + AVBipolar[u, v]
    }
  };


(* ::Subsection:: *)
(*Global styles for plots*)


textStyle = Style[#, 18] &;
textStyleBigger = Style[#, 20] &;


natStyle = Red;


flat0Style = Darker[Green];
sharp0Style = Magenta;


flatPiStyle = Blue;
sharpPiStyle = Orange;


nonStyle = Directive[Opacity[0.7], LightGray];
termStyle = Pink;
unphysStyle = Black;
inflStyle = Directive[Dashed, Green];


upperStyle = Blue;
lowerStyle = Red;
convexStyle = Black;


guideStyle = Dashed;
pointStyle = PointSize[Large];
glowStyle = Directive[Thick, Yellow, Opacity[0.7]];


uStyle = Blue;
vStyle = Red;


lighter = Directive[Opacity[0.5], #] &;


transStyle =
  Function[{a},
    If[a == aNat0 || a == aNatPi,
      Style[#, natStyle] &,
      Identity
    ]
  ];


(* ::Section:: *)
(*Bipolar coordinates*)


(* ::Subsection:: *)
(*Interactive visualiser*)


DynamicModule[
 {uInit, vInit, xyInit,
  vMax, xMax, yMax
 },
  uInit = 90 Degree;
  vInit = 3/10;
  xyInit = XYBipolar[uInit, vInit];
  vMax = 5;
  xMax = 4;
  yMax = 3;
  Manipulate[
    Show[
      EmptyFrame[{-xMax, xMax}, {-yMax, yMax},
        ImageSize -> 480,
        PlotLabel -> BoxedLabel @ Row[
          {
            Style[uIt == Row @ {u / Degree // N, Degree}, uStyle],
            Style[vIt == N[v], vStyle]
          }, ","
        ]
      ],
      (* All u-contours and v-contours *)
      uContours,
      vContours,
      (* Current u-contour and v-contour *)
      uContour[u, lighter @ uStyle],
      vContour[v, lighter @ vStyle],
      (* Current local orthonormal basis vectors *)
      uVector[u, v, uStyle],
      vVector[u, v, vStyle]
    ]
    (* NOTE: In Version 11, TrackingFunction won't work *)
    (* unless endpoint 2 Pi is numericised *)
    (* https://mathematica.stackexchange.com/a/101695 *)
  , {{u, uInit, uIt}, 0, 2 Pi,
      TrackingFunction -> Function[{uNew},
        u = uNew;
        xy = XYBipolar[u, v];
      ]
  }
  , {{v, vInit, vIt}, -vMax, vMax,
      TrackingFunction -> Function[{vNew},
        v = vNew;
        xy = XYBipolar[u, v];
      ]
  }
  , {{xy, xyInit}, Locator,
      TrackingFunction -> Function[{xyNew},
        xy = xyNew;
        u = UBipolar @@ xy // Mod[#, 2 Pi] &;
        v = VBipolar @@ xy;
      ]
  }]
]


(* ::Subsection:: *)
(*Animation (u varying, v fixed)*)


Module[
 {uValues, v,
  xMax, yMax
 },
  uValues = Subdivide[360, 36] Degree // Most;
  v = 0.3;
  xMax = 4;
  yMax = 3;
  Table[
    Show[
      EmptyFrame[{-xMax, xMax}, {-yMax, yMax},
        ImageSize -> 480,
        PlotLabel -> BoxedLabel @ Row[
          {
            Style[uIt == Row @ {u / Degree // N, Degree}, uStyle],
            Style[vIt == N[v], vStyle]
          }, ","
        ]
      ],
      (* All u-contours and v-contours *)
      uContours,
      vContours,
      (* Current u-contour and v-contour *)
      uContour[u, lighter @ uStyle],
      vContour[v, lighter @ vStyle],
      (* Current local orthonormal basis vectors *)
      uVector[u, v, uStyle],
      vVector[u, v, vStyle]
    ]
  , {u, uValues}]
] // Ex["bipolar-sweep-u.gif", gifOpts]


(* ::Subsection:: *)
(*Animation (u fixed, v varying)*)


Module[
 {u, vValues,
  xMax, yMax
 },
  u = 90 Degree;
  vValues = Subdivide[-4, 4, 20];
  xMax = 4;
  yMax = 3;
  Table[
    Show[
      EmptyFrame[{-xMax, xMax}, {-yMax, yMax},
        ImageSize -> 480,
        PlotLabel -> BoxedLabel @ Row[
          {
            Style[uIt == Row @ {u / Degree // N, Degree}, uStyle],
            Style[vIt == N[v], vStyle]
          }, ","
        ]
      ],
      (* All u-contours and v-contours *)
      uContours,
      vContours,
      (* Current u-contour and v-contour *)
      uContour[u, lighter @ uStyle],
      vContour[v, lighter @ vStyle],
      (* Current local orthonormal basis vectors *)
      uVector[u, v, uStyle],
      vVector[u, v, vStyle]
    ]
  , {v, vValues}]
] // Ex["bipolar-sweep-v.gif", gifOpts]


(* ::Section:: *)
(*Critical terminal points*)


Module[
 {aCold, aWarm, aHot,
  vMin, vMax,
  regList,
  aAss,
  a, nString,
  pointAlongCurve, lineToAxis, textBelowAxis
 },
  (* Non-borderline representative values of A *)
  aCold = Ceiling[1.1 aNat0, 1/2];
  aWarm = Way[aNatPi, aNat0];
  aHot = Floor[0.95 aNatPi, 1/2];
  (* Horizontal plot range *)
  vMin = Floor[vFlat0[aHot], 1/2];
  vMax = Ceiling[vSharp0[aHot], 1/2];
  (* List of regimes *)
  regList = {"hot", "warm_hot", "warm", "cold_warm", "cold"};
  (* Associations for various quantities *)
  aAss = AssociationThread[regList -> {aHot, aNatPi, aWarm, aNat0, aCold}];
  (* Plots for all five regimes *)
  Table[
    a = aAss[reg];
    nString = ToString @ First @ FirstPosition[regList, reg];
    pointAlongCurve = Point @ {#, #^4 / a} &;
    lineToAxis = Line @ {{#, #^4 / a}, {#, 0}} &;
    textBelowAxis = Function[{v, sub},
      Text[
        Subscript[
          vIt,
          sub // PrettyString[
            "Sharp" -> "\[Sharp]",
            "Nat" -> "\[Natural]",
            "Flat" -> "\[Flat]",
            "Pi" -> "\[Pi]"
          ]
        ] // textStyle,
        {v, 0},
        {0, 2.5}
      ]
    ];
    Show[
      (* Functions of v *)
      Plot[{Cosh[v] + 1, Cosh[v] - 1, v^4 / a}, {v, vMin, vMax},
        AxesLabel -> {vIt, Null},
        ImageSize -> 640,
        PlotLabel -> BoxedLabel[
          Row[
            {aIt == N[a], reg // PrettyString["_" -> "-to-"] // ""},
            "  "
          ] // transStyle[a]
        ],
        PlotLegends -> {
          Row[{Cosh[vIt] + 1, uIt == Pi // ""}, "  "],
          Row[{Cosh[vIt] - 1, uIt == 0 // ""}, "  "],
          Row[{vIt^4, aIt}, "/"]
        },
        PlotRange -> All,
        PlotRangeClipping -> False,
        PlotStyle -> {lighter[Blue], lighter[Red], Black},
        PlotOptions[Axes] // Evaluate
      ],
      (* Critical terminal points along u == 0 *)
      Which[
        (* Two distinct terminal points *)
        a < aNat0,
        {
          (* v_\[Flat]0 *)
          Graphics @ {Directive[flat0Style, pointStyle],
            pointAlongCurve @ vFlat0[a]
          },
          Graphics @ {Directive[flat0Style, guideStyle],
            lineToAxis @ vFlat0[a]
          },
          Graphics @ {flat0Style,
            textBelowAxis[vFlat0[a], "Flat0"]
          },
          (* v_\[Sharp]0 *)
          Graphics @ {Directive[sharp0Style, pointStyle],
            pointAlongCurve @ vSharp0[a]
          },
          Graphics @ {Directive[sharp0Style, guideStyle],
            lineToAxis @ vSharp0[a]
          },
          Graphics @ {sharp0Style,
            textBelowAxis[vSharp0[a], "Sharp0"]
          }
        },
        (* One terminal point *)
        a == aNat0,
        {
          (* v_\[Natural]0 *)
          Graphics @ {Directive[natStyle, pointStyle],
            pointAlongCurve @ vNat0
          },
          Graphics @ {Directive[natStyle, guideStyle],
            lineToAxis @ vNat0
          },
          Graphics @ {natStyle,
            textBelowAxis[vNat0, "Nat0"]
          }
        },
        (* Zero terminal points *)
        True,
        {}
      ],
      (* Critical terminal points along u == \[Pi] *)
      Which[
        (* Two distinct terminal points *)
        a < aNatPi,
        {
          (* v_\[Flat]\[Pi] *)
          Graphics @ {Directive[flatPiStyle, pointStyle],
            pointAlongCurve @ vFlatPi[a]
          },
          Graphics @ {Directive[flatPiStyle, guideStyle],
            lineToAxis @ vFlatPi[a]
          },
          Graphics @ {flatPiStyle,
            textBelowAxis[vFlatPi[a], "FlatPi"]
          },
          (* v_\[Sharp]\[Pi] *)
          Graphics @ {Directive[sharpPiStyle, pointStyle],
            pointAlongCurve @ vSharpPi[a]
          },
          Graphics @ {Directive[sharpPiStyle, guideStyle],
            lineToAxis @ vSharpPi[a]
          },
          Graphics @ {sharpPiStyle,
            textBelowAxis[vSharpPi[a], "SharpPi"]
          }
        },
        (* One terminal point *)
        a == aNatPi,
        {
          (* v_\[Natural]\[Pi] *)
          Graphics @ {Directive[natStyle, pointStyle],
            pointAlongCurve @ vNatPi
          },
          Graphics @ {Directive[natStyle, guideStyle],
            lineToAxis @ vNatPi
          },
          Graphics @ {natStyle,
            textBelowAxis[vNatPi, "NatPi"]
          }
        },
        (* Zero terminal points *)
        True,
        {
          Graphics @ {Transparent (* for consistent plot height *),
            textBelowAxis[vNatPi, "NatPi"]
          }
        }
      ]
    ] // Ex @ StringJoin["bipolar-critical-", nString, "-", reg, ".pdf"]
  , {reg, regList}]
]


(* ::Section:: *)
(*Viable domain animations*)


Module[
 {aMin, aMax, aStep,
  aValues,
  xMax, yMax,
  eps, xMaxUnphys, yMaxUnphys,
  xMinNon, xMaxNon, yMaxNon
 },
  (* Values of A *)
  aMin = 1;
  aMax = 10;
  aStep = 1/4;
  aValues = Range[aMin, aMax, aStep];
  aValues = Join[aValues, {aNat0, aNatPi}] // Sort[#, Less] &;
  (* Animation *)
  Table[
    (* Plot range *)
    xMax = 3.5;
    yMax = 2.5;
    (* Plot range for unphysical domain *)
    eps = 0.1;
    xMaxUnphys = xMax (1 + eps);
    yMaxUnphys = yMax (1 + eps);
    (* Plot range for non-viable domain *)
    (* (these bounds are merely guesses which work) *)
    xMinNon = If[a < aNatPi, vFlatPi[a], vNatPi] // 0.9 XBipolar[Pi, #] &;
    xMaxNon = If[a < aNat0, vFlat0[a], vNat0] // 1.01 XBipolar[0, #] &;
    yMaxNon = (xMaxNon - xMinNon) / 2;
    (* Plot *)
    Show[
      EmptyFrame[{-xMax, xMax}, {-yMax, yMax},
        PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
      ],
      (* Unphysical domain *)
      RegionPlot[x < 0,
        {x, -xMaxUnphys, eps}, {y, -yMaxUnphys, yMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ],
      (* All u-contours and v-contours *)
      uContours,
      vContours,
      (* Non-viable domain *)
      RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
        {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
        BoundaryStyle -> termStyle,
        PlotStyle -> nonStyle
      ]
    ]
  , {a, aValues}]
] // Ex["bipolar-viable-full.gif", gifOpts]


Module[
 {aMin, aMax, aStep,
  aValues,
  pointAlongAxis, textAboveBelowAxis,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon
 },
  (* Values of A *)
  aMin = 8;
  aMax = 10;
  aStep = 1/10;
  aValues = Range[aMin, aMax, aStep];
  aValues = Join[aValues, {aNat0, aNatPi}] // Sort[#, Less] &;
  (* Graphics abbreviations *)
  pointAlongAxis[u_, v_] := Point @ {XYBipolar[u, v]};
  textAboveBelowAxis[u_, v_, sub_, ysign_: 1] :=
    Text[
      Subscript[
        vIt,
        sub // PrettyString[
          "Sharp" -> "\[Sharp]",
          "Nat" -> "\[Natural]",
          "Flat" -> "\[Flat]",
          "Pi" -> "\[Pi]"
        ]
      ] // textStyle // BoxedLabel[#, Background -> White] &,
      XYBipolar[u, v],
      {0, -Sign[ysign] 1.7}
    ];
  (* Animation *)
  Table[
    (* Plot range *)
    xMin = vFlatPi[aMin] // 0.99 XBipolar[Pi, #] &;
    xMax = vFlat0[aMin] // 1.01 XBipolar[0, #] &;
    yMax = (xMax - xMin) / 2;
    (* Plot range for non-viable domain *)
    (* (these bounds are merely guesses which work) *)
    xMinNon = If[a < aNatPi, vFlatPi[a], vNatPi] // 0.9 XBipolar[Pi, #] &;
    xMaxNon = If[a < aNat0, vFlat0[a], vNat0] // 1.01 XBipolar[0, #] &;
    yMaxNon = (xMaxNon - xMinNon) / 2;
    (* Plot *)
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        ImageSize -> 480,
        PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
      ],
      (* All u-contours and v-contours *)
      uContours,
      vContours,
      (* Non-viable domain *)
      RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
        {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
        BoundaryStyle -> termStyle,
        PlotPoints -> 50,
        PlotStyle -> nonStyle
      ],
      (* Critical terminal points along u == 0 *)
      Which[
        (* Two distinct terminal points *)
        a < aNat0,
        {
          (* v_\[Flat]0 *)
          Graphics @ {Directive[flat0Style, pointStyle],
            pointAlongAxis[0, vFlat0[a]]
          },
          Graphics @ {flat0Style,
            textAboveBelowAxis[0, vFlat0[a], "Flat0", -1]
          },
          (* v_\[Sharp]0 *)
          Graphics @ {Directive[sharp0Style, pointStyle],
            pointAlongAxis[0, vSharp0[a]]
          },
          Graphics @ {sharp0Style,
            textAboveBelowAxis[0, vSharp0[a], "Sharp0"]
          }
        },
        (* One terminal point *)
        a == aNat0,
        {
          (* v_\[Natural]0 *)
          Graphics @ {Directive[natStyle, pointStyle],
            pointAlongAxis[0, vNat0]
          },
          Graphics @ {natStyle,
            textAboveBelowAxis[0, vNat0, "Nat0"]
          }
        },
        (* Zero terminal points *)
        True,
        {}
      ],
      (* Critical terminal points along u == \[Pi] *)
      Which[
        (* Two distinct terminal points *)
        a < aNatPi,
        {
          (* v_\[Flat]\[Pi] *)
          Graphics @ {Directive[flatPiStyle, pointStyle],
            pointAlongAxis[Pi, vFlatPi[a]]
          },
          Graphics @ {flatPiStyle,
            textAboveBelowAxis[Pi, vFlatPi[a], "FlatPi"]
          },
          (* v_\[Sharp]\[Pi] *)
          Graphics @ {Directive[sharpPiStyle, pointStyle],
            pointAlongAxis[Pi, vSharpPi[a]]
          },
          Graphics @ {sharpPiStyle,
            textAboveBelowAxis[Pi, vSharpPi[a], "SharpPi", -1]
          }
        },
        (* One terminal point *)
        a == aNatPi,
        {
          (* v_\[Natural]\[Pi] *)
          Graphics @ {Directive[natStyle, pointStyle],
            pointAlongAxis[Pi, vNatPi]
          },
          Graphics @ {natStyle,
            textAboveBelowAxis[Pi, vNatPi, "NatPi"]
          }
        },
        (* Zero terminal points *)
        True,
        {}
      ]
    ]
  , {a, aValues}]
] // Ex["bipolar-viable-zoom.gif", gifOpts]


(* ::Section:: *)
(*Traced boundary plots*)


(* ::Subsection:: *)
(*Hot regime A < A_\[Natural]\[Pi]*)


(* ::Subsubsection:: *)
(*Starting points*)


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon,
  idList
 },
  a = aHot;
  (* Plot range *)
  xMin = vFlatPi[a] // 0.99 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 1.07 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vFlatPi[a] // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Group names *)
  idList = {"outer", "inner", "axis", "hyperbolic"};
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 480,
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Starting points *)
    ListPlot[
      Table[
        XYBipolar @@@ startUVHot[id]
      , {id, idList}],
      LabelingFunction -> Function @ Placed[
        #2[[2]],
        If[#2[[1]] == 4, Left, Center]
      ],
      PlotLegends -> idList,
      PlotStyle -> Directive[Opacity[0.7], pointStyle]
    ]
  ]
] // Ex["bipolar-hot-traced-starting.pdf"]


(* ::Subsubsection:: *)
(*General*)


Module[
 {a,
  xMin, xMax, yMax,
  eps, xMaxUnphys, yMaxUnphys,
  xMinNon, xMaxNon, yMaxNon
 },
  a = aHot;
  (* Plot range *)
  xMin = vFlatPi[a] // 0 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 2 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for unphysical domain *)
  eps = 0.1;
  xMaxUnphys = xMax (1 + eps);
  yMaxUnphys = yMax (1 + eps);
  (* Plot range for non-viable domain *)
  xMinNon = vFlatPi[a] // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* Unphysical domain *)
    RegionPlot[x < 0,
      {x, -xMaxUnphys, eps}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> None,
      PlotStyle -> unphysStyle
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraHot[id]}]
    , {id, {"outer", "inner", "axis"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraHot[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-hot-traced-full.pdf"]


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon
 },
  a = aHot;
  (* Plot range *)
  xMin = vFlatPi[a] // 0.95 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 1.05 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vFlatPi[a] // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraHot[id]}]
    , {id, {"outer", "inner", "axis"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraHot[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-hot-traced-zoom_1.pdf"]


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon
 },
  a = aHot;
  (* Plot range *)
  xMin = vSharpPi[a] // 0.995 XBipolar[Pi, #] &;
  xMax = vSharp0[a] // 1.005 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vSharpPi[a] // 0.99 XBipolar[Pi, #] &;
  xMaxNon = vSharp0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraHot[id]}]
    , {id, {"inner"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraHot[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-hot-traced-zoom_2.pdf"]


(* ::Subsubsection:: *)
(*Convex without known solution*)


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon,
  v
 },
  a = aHot;
  (* Plot range *)
  xMin = vSharpPi[a] // 0.995 XBipolar[Pi, #] &;
  xMax = vSharp0[a] // 1.005 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vSharpPi[a] // 0.99 XBipolar[Pi, #] &;
  xMaxNon = vSharp0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Traced boundaries *)
  v = vTraHot["hyperbolic"] // First;
  v = Function[{u}, v[-Abs @ u] // Evaluate];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    ParametricPlot[
      XYBipolar[u, v[u]],
      {u, -Pi, Pi},
      PlotStyle -> convexStyle
    ]
  ]
] // Ex["bipolar-hot-traced-convex.pdf"]


(* ::Subsubsection:: *)
(*Convex with known solution*)


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon,
  v, vBath
 },
  a = aHot;
  (* Plot range *)
  xMin = vSharpPi[a] // 0.995 XBipolar[Pi, #] &;
  xMax = vSharp0[a] // 1.005 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vSharpPi[a] // 0.99 XBipolar[Pi, #] &;
  xMaxNon = vSharp0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Traced boundaries *)
  v = vTraHot["hyperbolic"] // First;
  v = Function[{u}, v[-Abs @ u] // Evaluate];
  (* Heat bath *)
  vBath = 1.2 vSharp0[a];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Known solution *)
    DensityPlot[VBipolar[x, y],
      {x, xMin, xMax}, {y, -yMax, yMax},
      ColorFunction -> "TemperatureMap",
      RegionFunction -> Function[{x, y},
        v @ UBipolar[x, y] < VBipolar[x, y] < vBath
          // Evaluate
      ]
    ],
    (* Traced boundaries *)
    ParametricPlot[
      XYBipolar[u, v[u]],
      {u, -Pi, Pi},
      PlotStyle -> convexStyle
    ]
  ]
] // Ex["bipolar-hot-traced-convex-with-known.png"]


(* ::Subsection:: *)
(*Warm-to-hot transition A = A_\[Natural]\[Pi]*)


(* ::Subsubsection:: *)
(*Starting points*)


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon,
  idList
 },
  a = aWarmHot;
  (* Plot range *)
  xMin = vFlatPi[a] // 0.99 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 1.07 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Group names *)
  idList = {"outer", "inner", "axis", "hyperbolic"};
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 480,
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 50,
      PlotStyle -> nonStyle
    ],
    (* Starting points *)
    ListPlot[
      Table[
        XYBipolar @@@ startUVWarmHot[id]
      , {id, idList}],
      LabelingFunction -> Function @ Placed[
        #2[[2]],
        If[#2[[1]] == 4, Left, Center]
      ],
      PlotLegends -> idList,
      PlotStyle -> Directive[Opacity[0.7], pointStyle]
    ]
  ]
] // Ex["bipolar-warm_hot-traced-starting.pdf"]


(* ::Subsubsection:: *)
(*General*)


Module[
 {a,
  xMin, xMax, yMax,
  eps, xMaxUnphys, yMaxUnphys,
  xMinNon, xMaxNon, yMaxNon
 },
  a = aWarmHot;
  (* Plot range *)
  xMin = vNatPi // 0 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 2 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for unphysical domain *)
  eps = 0.1;
  xMaxUnphys = xMax (1 + eps);
  yMaxUnphys = yMax (1 + eps);
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* Unphysical domain *)
    RegionPlot[x < 0,
      {x, -xMaxUnphys, eps}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> None,
      PlotStyle -> unphysStyle
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraWarmHot[id]}]
    , {id, {"outer", "inner", "axis"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraWarmHot[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-warm_hot-traced-full.pdf"]


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon
 },
  a = aWarmHot;
  (* Plot range *)
  xMin = vNatPi // 0.95 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 1.05 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 50,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraWarmHot[id]}]
    , {id, {"outer", "inner", "axis"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraWarmHot[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-warm_hot-traced-zoom_1.pdf"]


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon
 },
  a = aWarmHot;
  (* Plot range *)
  xMin = vNatPi // 0.995 XBipolar[Pi, #] &;
  xMax = vSharp0[a] // 1.005 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.99 XBipolar[Pi, #] &;
  xMaxNon = vSharp0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 70,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraWarmHot[id]}]
    , {id, {"outer", "inner"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraWarmHot[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-warm_hot-traced-zoom_2.pdf"]


(* ::Subsubsection:: *)
(*Convex without known solution*)


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon,
  v
 },
  a = aWarmHot;
  (* Plot range *)
  xMin = vNatPi // 0.995 XBipolar[Pi, #] &;
  xMax = vSharp0[a] // 1.005 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.99 XBipolar[Pi, #] &;
  xMaxNon = vSharp0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Traced boundaries *)
  v = vTraWarmHot["hyperbolic"] // First;
  v = Function[{u}, v[-Abs @ u] // Evaluate];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 70,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    ParametricPlot[
      XYBipolar[u, v[u]],
      {u, -Pi, Pi},
      PlotStyle -> convexStyle
    ]
  ]
] // Ex["bipolar-warm_hot-traced-convex.pdf"]


(* ::Subsubsection:: *)
(*Convex with known solution*)


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon,
  v, vBath
 },
  a = aWarmHot;
  (* Plot range *)
  xMin = vNatPi // 0.995 XBipolar[Pi, #] &;
  xMax = vSharp0[a] // 1.005 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.99 XBipolar[Pi, #] &;
  xMaxNon = vSharp0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Traced boundaries *)
  v = vTraWarmHot["hyperbolic"] // First;
  v = Function[{u}, v[-Abs @ u] // Evaluate];
  (* Heat bath *)
  vBath = 1.2 vSharp0[a];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 70,
      PlotStyle -> nonStyle
    ],
    (* Known solution *)
    DensityPlot[VBipolar[x, y],
      {x, xMin, xMax}, {y, -yMax, yMax},
      ColorFunction -> "TemperatureMap",
      RegionFunction -> Function[{x, y},
        v @ UBipolar[x, y] < VBipolar[x, y] < vBath
          // Evaluate
      ]
    ],
    (* Traced boundaries *)
    ParametricPlot[
      XYBipolar[u, v[u]],
      {u, -Pi, Pi},
      PlotStyle -> convexStyle
    ]
  ]
] // Ex["bipolar-warm_hot-traced-convex-with-known.png"]


(* ::Subsection:: *)
(*Warm regime A_\[Natural]\[Pi] < A < A_\[Natural]0*)


(* ::Subsubsection:: *)
(*Starting points*)


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon,
  idList
 },
  a = aWarm;
  (* Plot range *)
  xMin = vNatPi // 0.99 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 1.07 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Group names *)
  idList = {"terminal", "axis", "hyperbolic"};
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 480,
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 50,
      PlotStyle -> nonStyle
    ],
    (* Starting points *)
    ListPlot[
      Table[
        XYBipolar @@@ startUVWarm[id]
      , {id, idList}],
      LabelingFunction -> Function @ Placed[
        #2[[2]],
        If[#2[[1]] == 3, Left, Center]
      ],
      PlotLegends -> idList,
      PlotStyle -> Directive[Opacity[0.7], pointStyle]
    ]
  ]
] // Ex["bipolar-warm-traced-starting.pdf"]


(* ::Subsubsection:: *)
(*General*)


Module[
 {a,
  xMin, xMax, yMax,
  eps, xMaxUnphys, yMaxUnphys,
  xMinNon, xMaxNon, yMaxNon
 },
  a = aWarm;
  (* Plot range *)
  xMin = vNatPi // 0 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 2 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for unphysical domain *)
  eps = 0.1;
  xMaxUnphys = xMax (1 + eps);
  yMaxUnphys = yMax (1 + eps);
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* Unphysical domain *)
    RegionPlot[x < 0,
      {x, -xMaxUnphys, eps}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> None,
      PlotStyle -> unphysStyle
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraWarm[id]}]
    , {id, {"terminal", "axis"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraWarm[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-warm-traced-full.pdf"]


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon
 },
  a = aWarm;
  (* Plot range *)
  xMin = vNatPi // 0.97 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 1.03 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 50,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraWarm[id]}]
    , {id, {"terminal", "axis"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraWarm[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-warm-traced-zoom.pdf"]


(* ::Subsubsection:: *)
(*Convex without known solution*)


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon,
  v
 },
  a = aWarm;
  (* Plot range *)
  xMin = vNatPi // 0.97 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vSharp0[a] // 1.1 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Traced boundaries *)
  v = vTraWarm["hyperbolic"] // First;
  v = Function[{u}, v[-Abs @ u] // Evaluate];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 50,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries *)
    ParametricPlot[
      XYBipolar[u, v[u]],
      {u, -Pi, Pi},
      PlotStyle -> convexStyle
    ]
  ]
] // Ex["bipolar-warm-traced-convex.pdf"]


(* ::Subsubsection:: *)
(*Convex with known solution*)


Module[
 {a,
  xMin, xMax, yMax,
  xMinNon, xMaxNon, yMaxNon,
  v, vBath
 },
  a = aWarm;
  (* Plot range *)
  xMin = vNatPi // 0.97 XBipolar[Pi, #] &;
  xMax = vFlat0[a] // 1.01 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for non-viable domain *)
  xMinNon = vNatPi // 0.9 XBipolar[Pi, #] &;
  xMaxNon = vSharp0[a] // 1.1 XBipolar[0, #] &;
  yMaxNon = (xMaxNon - xMinNon) / 2;
  (* Traced boundaries *)
  v = vTraWarm["hyperbolic"] // First;
  v = Function[{u}, v[-Abs @ u] // Evaluate];
  (* Heat bath *)
  vBath = 1.2 vSharp0[a];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Non-viable domain *)
    RegionPlot[vi[a] @@ UVBipolar[x, y] < 0,
      {x, xMinNon, xMaxNon}, {y, -yMaxNon, yMaxNon},
      BoundaryStyle -> termStyle,
      PlotPoints -> 50,
      PlotStyle -> nonStyle
    ],
    (* Known solution *)
    DensityPlot[VBipolar[x, y],
      {x, xMin, xMax}, {y, -yMax, yMax},
      ColorFunction -> "TemperatureMap",
      RegionFunction -> Function[{x, y},
        v @ UBipolar[x, y] < VBipolar[x, y] < vBath
          // Evaluate
      ]
    ],
    (* Traced boundaries *)
    ParametricPlot[
      XYBipolar[u, v[u]],
      {u, -Pi, Pi},
      PlotStyle -> convexStyle
    ]
  ]
] // Ex["bipolar-warm-traced-convex-with-known.png"]


(* ::Subsection:: *)
(*Cold-to-warm transition A = A_\[Natural]0*)


(* ::Subsubsection:: *)
(*Starting points*)


Module[
 {a,
  xMin, xMax, yMax,
  idList
 },
  a = aColdWarm;
  (* Plot range *)
  xMin = 0.99;
  xMax = vNat0 // 2.01 XBipolar[0, #] - 1 &;
  yMax = (xMax - xMin) / 2;
  (* Group names *)
  idList = {"axis", "hyperbolic"};
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 480,
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Starting points *)
    ListPlot[
      Table[
        XYBipolar @@@ startUVColdWarm[id]
      , {id, idList}],
      LabelingFunction -> Function @ Placed[
        #2[[2]],
        If[#2[[1]] == 2, Right, Center]
      ],
      PlotLegends -> idList,
      PlotStyle -> Directive[Opacity[0.7], pointStyle]
    ]
  ]
] // Ex["bipolar-cold_warm-traced-starting.pdf"]


(* ::Subsubsection:: *)
(*General*)


Module[
 {a,
  xMin, xMax, yMax,
  eps, xMaxUnphys, yMaxUnphys
 },
  a = aColdWarm;
  (* Plot range *)
  xMin = vNatPi // 0 XBipolar[Pi, #] &;
  xMax = vNat0 // 2 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for unphysical domain *)
  eps = 0.1;
  xMaxUnphys = xMax (1 + eps);
  yMaxUnphys = yMax (1 + eps);
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* Unphysical domain *)
    RegionPlot[x < 0,
      {x, -xMaxUnphys, eps}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> None,
      PlotStyle -> unphysStyle
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraColdWarm[id]}]
    , {id, {"axis"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraColdWarm[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-cold_warm-traced-full.pdf"]


Module[
 {a,
  xMin, xMax, yMax
 },
  a = aColdWarm;
  (* Plot range *)
  xMin = 6/8;
  xMax = 9/8;
  yMax = (xMax - xMin) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraColdWarm[id]}]
    , {id, {"axis"}}],
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> glowStyle
        ]
      , {v, vTraColdWarm[id]}]
    , {id, {"hyperbolic"}}]
  ]
] // Ex["bipolar-cold_warm-traced-zoom.pdf"]


(* ::Subsubsection:: *)
(*Pseudo-convex without known solution*)


Module[
 {a,
  xMin, xMax, yMax,
  v
 },
  a = aColdWarm;
  (* Plot range *)
  xMin = 6/8;
  xMax = 17/16;
  yMax = (xMax - xMin) / 2;
  (* Traced boundaries *)
  v = vTraColdWarm["hyperbolic"] // First;
  v = Function[{u}, v[-Abs @ u] // Evaluate];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Traced boundaries *)
    ParametricPlot[
      XYBipolar[u, v[u]],
      {u, -Pi, Pi},
      PlotStyle -> convexStyle
    ]
  ]
] // Ex["bipolar-cold_warm-traced-pseudo_convex.pdf"]


Module[
 {a,
  xMin, xMax, yMax,
  v, uStart,
  reflect, tanVec,
  uInfl, vInfl, xyInfl,
  tanVecInfl, tanLineInfl,
  xyCorn, xyInt,
  mark,
  newLine = "\n",
  concaveStyle = Red,
  inflDotStyle = Directive[Yellow, Opacity[0.5], pointStyle],
  inflMarkStyle = Directive[Thin, Darker[Green]],
  tanStyle = Directive[Blue, Thin],
  cornDotStyle = Directive[Black, Opacity[0.7], pointStyle],
  cornMarkStyle = Directive[Red, Thin],
  intDotStyle = Directive[Cyan, Opacity[0.5], pointStyle],
  intMarkStyle = Directive[Blue, Opacity[0.7], Thin]
 },
  a = aColdWarm;
  (* Plot range *)
  xMin = 6/8;
  xMax = 17/16;
  yMax = (xMax - xMin) / 2;
  (* Traced boundaries (NOTE: uStart is NEGATIVE) *)
  v = vTraColdWarm["hyperbolic"] // First;
  uStart = DomainStart[v];
  v = Function[{u}, v[-Abs @ u] // Evaluate];
  (* Reflection across x-axis *)
  reflect[{x_, y_}] := {x, -y};
  (* Cartesian components of tangent vector to traced boundary *)
  tanVec[u_, v_] :=
    HBipolar[u, v] Plus[
      AUBipolar[u, v],
      vTraDer[a][u, v] AVBipolar[u, v]
    ] // Evaluate;
  (* Points of inflection (NOTE: uInfl is NEGATIVE) *)
  uInfl = SeekRoot[curTra[a][#, v[#]] &, {-Pi, 0}];
  vInfl = v[uInfl];
  xyInfl = XYBipolar[uInfl, vInfl];
  tanVecInfl = tanVec[uInfl, vInfl];
  tanLineInfl = InfiniteLine[xyInfl, tanVecInfl];
  (* Corner formed by concave portions *)
  xyCorn = XYBipolar[-Pi, v[-Pi]];
  (* Intersection of tangent lines *)
  xyInt =
    With[{x = \[FormalX], y = \[FormalY]},
      {x, y} /. First @ FindInstance[
        And[
          Element[{x, y}, tanLineInfl],
          Element[{x, y}, reflect /@ tanLineInfl]
        ]
      , {x, y}]
    ];
  (* Marker *)
  mark[{x_, y_}, rot_, len_: 0.01] := (
    Line @ {{0, -len}, {0, len}}
      // Rotate[#, rot] &
      // Translate[#, {x, y}] &
  );
  (* Plots *)
  {
    (* Traced boundary convexity as a function of u *)
    Show[
      Plot[
        curTra[a][u, v[u]] {
          If[u < uInfl, 1, Indeterminate],
          If[u > uInfl, 1, Indeterminate]
        } // Evaluate, {u, -Pi, 0},
        AxesLabel -> {uIt, "Cur."},
        PlotRange -> Full,
        PlotStyle -> {concaveStyle, convexStyle},
        PlotOptions[Axes] // Evaluate
      ],
      (* Inflection *)
      Graphics @ {inflDotStyle,
        Point @ {uInfl, 0}
      },
      Graphics @ Text[
        Subscript[uIt, "i"] == uInfl // textStyle,
        {uInfl, 0},
        {-1.1, -1.5}
      ]
    ] // Ex["bipolar-cold_warm-traced-pseudo_convex-curvature.pdf"],
    (* Detailed plot demonstrating inflection and concavity *)
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax},
        PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
      ],
      (* All u-contours and v-contours *)
      uContours,
      vContours,
      (* Traced boundaries (convex portion) *)
      ParametricPlot[
        XYBipolar[u, v[u]],
        {u, uInfl, -uInfl},
        PlotStyle -> convexStyle
      ],
      (* Traced boundaries (concave portions) *)
      ParametricPlot[
        XYBipolar[u, v[u]],
        {u, uStart, uInfl},
        PlotStyle -> concaveStyle
      ],
      ParametricPlot[
        XYBipolar[u, v[u]],
        {u, -uInfl, -uStart},
        PlotStyle -> concaveStyle
      ],
      (* Points of inflection *)
      Graphics @ {inflDotStyle,
        Point @ {xyInfl, xyInfl // reflect}
      },
      (* Mark at points of inflection *)
      Graphics @ {inflMarkStyle,
        mark[xyInfl, {{1, 0}, tanVecInfl}],
        mark[xyInfl // reflect, {{1, 0}, tanVecInfl // reflect}]
      },
      (* Tangent line at points of inflection *)
      Graphics @ {tanStyle,
        tanLineInfl,
        reflect /@ tanLineInfl
      },
      (* Corner formed by concave portions *)
      Graphics @ {cornDotStyle,
        Point @ {xyCorn}
      },
      (* Marker at corner *)
      Graphics @ {cornMarkStyle,
        mark[xyCorn, 0]
      },
      (* Intersection of tangent lines *)
      Graphics @ {intDotStyle,
        Point @ {xyInt}
      },
      (* Mark at intersection of tangent lines *)
      Graphics @ {intMarkStyle,
        mark[xyInt, 0]
      },
      (* Explanatory text *)
      Graphics @ {
        Text[
          StringJoin[
            "Black: convex portion of boundary", newLine,
            "Red: concave portions of boundary", newLine,
            "Green marks (yellow dots): points of inflection", newLine,
            "Blue: tangent line at points of inflection", newLine,
            "Red mark (black dot): corner formed by concave portions", newLine,
            "Blue mark (cyan dot): intersection of tangent lines", newLine,
            newLine,
            "This PDF consists of vector graphics.", newLine,
            "Zoom in to see the small gap between the red and blue marks."
          ] // Style[#, 7, Background -> White] &,
          {Way[xMin, xMax], Way[0, yMax, 7/10]}
        ]
      }
    ] // Ex["bipolar-cold_warm-traced-pseudo_convex-detailed.pdf"]
  }
]


(* ::Subsubsection:: *)
(*Pseudo-convex with known solution*)


Module[
 {a,
  xMin, xMax, yMax,
  v, vBath
 },
  a = aColdWarm;
  (* Plot range *)
  xMin = 6/8;
  xMax = 17/16;
  yMax = (xMax - xMin) / 2;
  (* Traced boundaries *)
  v = vTraColdWarm["hyperbolic"] // First;
  v = Function[{u}, v[-Abs @ u] // Evaluate];
  (* Heat bath *)
  vBath = 1.2 vSharp0[a];
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Known solution *)
    DensityPlot[VBipolar[x, y],
      {x, xMin, xMax}, {y, -yMax, yMax},
      ColorFunction -> "TemperatureMap",
      RegionFunction -> Function[{x, y},
        v @ UBipolar[x, y] < VBipolar[x, y] < vBath
          // Evaluate
      ]
    ],
    (* Traced boundaries *)
    ParametricPlot[
      XYBipolar[u, v[u]],
      {u, -Pi, Pi},
      PlotStyle -> convexStyle
    ]
  ]
] // Ex["bipolar-cold_warm-traced-pseudo_convex-with-known.png"]


(* ::Subsection:: *)
(*Cold regime A > A_\[Natural]0*)


(* ::Subsubsection:: *)
(*Starting points*)


Module[
 {a,
  xMin, xMax, yMax,
  idList
 },
  a = aCold;
  (* Plot range *)
  xMin = 1;
  xMax = 2;
  yMax = (xMax - xMin) / 2;
  (* Group names *)
  idList = {"axis"};
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      ImageSize -> 480,
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Starting points *)
    ListPlot[
      Table[
        XYBipolar @@@ startUVCold[id]
      , {id, idList}],
      LabelingFunction -> Function @ Placed[
        #2[[2]],
        If[#2[[1]] == 2, Right, Center]
      ],
      PlotLegends -> idList,
      PlotStyle -> Directive[Opacity[0.7], pointStyle]
    ]
  ]
] // Ex["bipolar-cold-traced-starting.pdf"]


(* ::Subsubsection:: *)
(*General*)


Module[
 {a,
  xMin, xMax, yMax,
  eps, xMaxUnphys, yMaxUnphys
 },
  a = aCold;
  (* Plot range *)
  xMin = vNatPi // 0 XBipolar[Pi, #] &;
  xMax = vNat0 // 2 XBipolar[0, #] &;
  yMax = (xMax - xMin) / 2;
  (* Plot range for unphysical domain *)
  eps = 0.1;
  xMaxUnphys = xMax (1 + eps);
  yMaxUnphys = yMax (1 + eps);
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* Unphysical domain *)
    RegionPlot[x < 0,
      {x, -xMaxUnphys, eps}, {y, -yMaxUnphys, yMaxUnphys},
      BoundaryStyle -> None,
      PlotStyle -> unphysStyle
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraCold[id]}]
    , {id, {"axis"}}]
  ]
] // Ex["bipolar-cold-traced-full.pdf"]


Module[
 {a,
  xMin, xMax, yMax
 },
  a = aCold;
  (* Plot range *)
  xMin = 0.95;
  xMax = 1.25;
  yMax = (xMax - xMin) / 2;
  (* Plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax},
      PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
    ],
    (* All u-contours and v-contours *)
    uContours,
    vContours,
    (* Traced boundaries *)
    Table[
      Table[
        ParametricPlot[
          {XYBipolar[u, v[u]], XYBipolar[-u, v[u]]},
          {u, DomainStart[v], DomainEnd[v]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {v, vTraCold[id]}]
    , {id, {"axis"}}]
  ]
] // Ex["bipolar-cold-traced-zoom.pdf"]


(* ::Section:: *)
(*Traced boundary algebra*)


(* ::Text:: *)
(*See Section 4.4 (Page r3-10).*)


(* ::Subsection:: *)
(*Check derivatives*)


(* ::Text:: *)
(*See (r3.56) to (r3.59)(Page r3-10).*)


With[{u = \[FormalU], v = \[FormalV], a = \[FormalCapitalA]},
  Module[
   {au, av, d, ss, cc, h, dd, ee
   },
    (* Local orthogonal basis *)
    au = AUBipolar[u, v];
    av = AVBipolar[u, v];
    (* Abbreviation for u-derivative *)
    d = Dt[#, u, Constants -> a] &;
    (* Abbreviations *)
    ss = Sin[u] Sinh[v];
    cc = Cos[u] Cosh[v] - 1;
    h = HBipolar[u, v];
    dd = Sinh[v] - d[v] Sin[u];
    ee = d[v] Sinh[v] + Sin[u];
    (* Table of results *)
    {
      {
        "(h S)' == h^2 C D",
        "(r3.56)",
        d[h ss] == h^2 cc dd // FullSimplify
      },
      {
        "(h C)' == -h^2 S D",
        "(r3.57)",
        d[h cc] == -h^2 ss dd // FullSimplify
      },
      {
        "(a_u)' == h D a_v",
        "(r3.58)",
        d[au] == h dd av // FullSimplify
      },
      {
        "(a_v)' == -h D a_u",
        "(r3.59)",
        d[av] == -h dd au // FullSimplify
      }
    } // TableForm
  ]
] // Ex["bipolar-check-derivatives.pdf"]


(* ::Subsection:: *)
(*Check curvature*)


(* ::Text:: *)
(*See (r3.62) and (r3.63) (Page r3-10).*)


With[{u = \[FormalU], v = \[FormalV], a = \[FormalCapitalA]},
  Module[
   {au, av, d, ss, cc, h, dd, ee,
    dSimp, vel, acc, comp, cur,
    newLine = "\n  "
   },
    (* Local orthogonal basis *)
    {au, av};
    (* Abbreviation for u-derivative *)
    d = Dt[#, u, Constants -> a] &;
    (* Abbreviations *)
    ss = Sin[u] Sinh[v];
    cc = Cos[u] Cosh[v] - 1;
    h = HBipolar[u, v];
    dd = Sinh[v] - d[v] Sin[u];
    ee = d[v] Sinh[v] + Sin[u];
    (* Simplification rules *)
    dSimp = {
      d[au] -> h dd av,
      d[av] -> -h dd au
    };
    (* Velocity and acceleration *)
    vel = h au + h d[v] av;
    acc = d[vel] /. dSimp;
    (* Components for local (u, v, z) basis *)
    comp = {
      au -> {1, 0, 0},
      av -> {0, 1, 0}
    };
    (* Curvature *)
    cur = Last @ Cross[
      vel /. comp,
      acc /. comp
    ];
    (* Table of results *)
    {
      {
        StringJoin[
          "Acc. ==", newLine,
          "- h^2 (D v' + E) a_u", newLine,
          "+ h (v'' - h E v' + h D) a_v"
        ],
        "(r3.62)",
        acc ==
          - h^2 (dd d[v] + ee) au
          + h (d @ d[v] - h ee d[v] + h dd) av
          /. comp
          // FullSimplify
      },
      {
        "Cur. == h^3 (D (1 + v'^2) + v'' / h)",
        "(r3.63)",
        cur == h^3 (dd (1 + d[v]^2) + d @ d[v] / h)
          /. comp 
          // FullSimplify
      }
    } // TableForm
  ]
] // Ex["bipolar-check-curvature.pdf"]


(* ::Subsection:: *)
(*Curvature at u = -\[Pi]*)


(* ::Text:: *)
(*See (r3.64) to (r3.66) (Page r3-10).*)


With[{v = \[FormalV], a = \[FormalCapitalA]},
  {
    {
      "Cur. (u = -\[Pi])",
      "(r3.65)",
      curTra[a][-Pi, v] // FullSimplify
    },
    {
      "v_i",
      "(r3.66)",
      vInfl // N
    },
    {
      "",
      "(check)",
      curTra[a][-Pi, vInfl] == 0 // FullSimplify
    }
  } // TableForm
] // Ex["bipolar-curvature-corner.pdf"]


(* ::Subsection:: *)
(*Hyperbolic coordinate at u = -\[Pi] for candidate boundaries*)


Module[
 {aMin, aMax, num, aValues,
  dest, aVTable,
  inflDotStyle = Directive[Red, Opacity[0.7], pointStyle]
 },
  (* Values of dimensionless group A *)
  (* (argument to vTraCandCorn must be EXACT) *)
  aMin = aNat0 - 4/1000;
  aMax = aNat0;
  num = 100;
  aValues = Subdivide[aMin, aMax, num];
  (* Compute table of (A, v(-\[Pi])) values *)
  (* (This is really slow (~25 sec), so compute once and store.) *)
  (* (Delete the file manually to compute from scratch.) *)
  dest = "bipolar-candidate-v-corner.txt";
  If[FileExistsQ[dest],
    (* If already stored, import *)
    aVTable = Import[dest] // Uncompress,
    (* Otherwise compute and export for next time *)
    aVTable = Table[
      {a // N, vTraCandCorn[a]}
    , {a, aValues}];
    aVTable // Compress // Ex[dest]
  ];
  (* Plot *)
  Show[
    (* A vs v(-\[Pi]) *)
    ListPlot[
      aVTable,
      AxesLabel -> {aIt, vIt[-Pi]},
      Joined -> True,
      PlotRange -> All,
      PlotStyle -> Black,
      PlotOptions[Axes] // Evaluate
    ],
    (* v_i *)
    Plot[vInfl, {a, aMin, aMax},
      PlotStyle -> inflStyle
    ],
    (* A_i *)
    Graphics @ {inflDotStyle,
      Point @ {aInfl, vInfl}
    },
    Graphics @ Text[
      Column @ {
        Subscript[aIt, "i"] == N[aInfl],
        Subscript[vIt, "i"] == N[vInfl]
      } // textStyle,
      {aInfl, vInfl},
      {-1.2, -1.2}
    ]
  ]
] // Ex["bipolar-candidate-v-corner.pdf"]


(* ::Section:: *)
(*Physical temperature range algebra*)


(* ::Text:: *)
(*See Section 5.1 (Page r3-14).*)


(* ::Subsection:: *)
(*\[Omega] properties*)


With[{v = \[FormalV]},
  Module[{omegaAsyInfinity},
    omegaAsyInfinity[v_] := 1 / v - 2 Exp[-v] / v;
    {
      {"omega", omega[v]},
      {"Small v series", omega[v] + O[v]^3},
      {"Large v asymptotic", omegaAsyInfinity[v]},
      {
        "(check)",
        And @@ (
          Limit[#, v -> Infinity] == 1 & /@
          {
            (1 / v) / omega[v],
            (omega[v] - (1 / v)) / (- 2 Exp[-v] / v)
          }
        )
      },
      {"domega/dv", omega'[v] // FullSimplify}
    } /. {v -> Subscript[v, "sharp0"]}
      // PrettyString["omega" -> "\[Omega]", "v" -> ToString[v] <> "_sharp0"]
      // PrettyString["sharp" -> "\[Sharp]"]
      // TableForm
  ]
] // Ex["bipolar-omega-properties.pdf"]


(* ::Subsection:: *)
(*\[Omega] plot*)


Module[{omegaAsyOrigin, omegaAsyInfinity, vItSharp0},
  omegaAsyOrigin[v_] := 1 / 2 - v^2 / 24;
  omegaAsyInfinity[v_] := 1 / v - 2 Exp[-v] / v;
  vItSharp0 = Subscript[vIt, "sharp0"];
  Show[
    Plot[
      {
        omega[v],
        Piecewise @ {
          {omegaAsyOrigin[v], v < 3.4},
          {Indeterminate, True}
        },
        Piecewise @ {
          {omegaAsyInfinity[v], v > 0.7},
          {Indeterminate, True}
        }
      }, {v, 0, 10},
      AxesLabel -> {vItSharp0, "omega"},
      PlotLegends -> {
        "omega"[vItSharp0],
        omegaAsyOrigin[vItSharp0],
        omegaAsyInfinity[vItSharp0]
      },
      PlotRange -> All,
      PlotRangeClipping -> False,
      PlotStyle -> {
        Blue,
        Directive[guideStyle, Darker[Green]],
        Directive[guideStyle, Magenta]
      },
      PlotOptions[Axes] // Evaluate
    ],
    (* {vNat0, 1/4} *)
    Graphics @ {Directive[guideStyle, natStyle],
      Line @ {{0, 1/4}, {vNat0, 1/4}, {vNat0, 0}}
    },
    Graphics @ {natStyle,
      Text[Row[{1, 4}, "/"] // textStyle, {0, 1/4}, {1.3, 0}],
      Text[Subscript[vIt, "nat0"] // textStyle, {vNat0, 0}, {0, 2}]
    }
  ]
    // PrettyString[
      "omega" -> "\[Omega]",
      "nat" -> "\[Natural]",
      "sharp" -> "\[Sharp]"
    ]
] // Ex["bipolar-omega.pdf"]


(* ::Section:: *)
(*Physical asymmetry*)


(* ::Subsection:: *)
(*Check asymmetry*)


(* ::Text:: *)
(*See (r3.92) (Page r3-16).*)


With[{v = \[FormalV]},
  Module[{vSharp0, vCorn},
    vSharp0 = Subscript[v, "sharp0"];
    vCorn = v[-Pi];
    {
      {"Asymmetry", asymm[vSharp0, vCorn]},
      {
        "(r3.92)",
        Equal[
          asymm[vSharp0, vCorn],
          Divide[
            1 - Sinh[vCorn] / (Cosh[vCorn] + 1),
            Sinh[vSharp0] / (Cosh[vSharp0] - 1) - 1
          ]
        ] // FullSimplify
      }
    }
      // PrettyString["sharp" -> "\[Sharp]"]
      // TableForm
  ]
] // Ex["bipolar-check-asymmetry.pdf"]


(* ::Subsection:: *)
(*Plot asymmetry*)


Module[
 {aMin, aMax, num, aValues,
  dest, aAsymmTable, asymmNatPi, asymmNat0
 },
  (* Values of dimensionless group A *)
  (* (argument to vTraCandCorn must be EXACT) *)
  aMin = 8;
  aMax = aNat0;
  num = 100;
  aValues = Subdivide[aMin, aMax, num];
  (* Compute table of (A, asymmetry) values *)
  (* (This is really slow (~35 sec), so compute once and store.) *)
  (* (Delete the file manually to compute from scratch.) *)
  dest = "bipolar-candidate-asymmetry.txt";
  If[FileExistsQ[dest],
    (* If already stored, import *)
    {aAsymmTable, asymmNatPi, asymmNat0} = Import[dest] // Uncompress,
    (* Otherwise compute and export for next time *)
    aAsymmTable = Table[
      {a // N, asymm[vSharp0[a], vTraCandCorn[a]]}
    , {a, aValues}];
    asymmNatPi = asymm[vSharp0[aNatPi], vTraCandCorn[aNatPi]];
    asymmNat0 = Last[aAsymmTable][[2]];
    {aAsymmTable, asymmNatPi, asymmNat0} // Compress // Ex[dest]
  ];
  (* Plot *)
  Show[
    (* A vs asymmetry *)
    ListPlot[
      aAsymmTable,
      AxesLabel -> {aIt, "Asymm."},
      Joined -> True,
      PlotRange -> All,
      PlotStyle -> Black,
      PlotOptions[Axes] // Evaluate
    ],
    (* aNatPi *)
    Graphics @ {Directive[natStyle, guideStyle],
      Line @ {{aMin, asymmNatPi}, {aNatPi, asymmNatPi}, {aNatPi, 0}}
    },
    Graphics @ {natStyle,
      Text[N[asymmNatPi], {aMin, asymmNatPi}, {1.2, 0}],
      Text[Subscript[aIt, "natpi"] // textStyle, {aNatPi, 0}, {0, 2.3}]
    },
    (* aNat0 *)
    Graphics @ {Directive[natStyle, guideStyle],
      Line @ {{aMin, asymmNat0}, {aNat0, asymmNat0}, {aNat0, 0}}
    },
    Graphics @ {natStyle,
      Text[N[asymmNat0], {aMin, asymmNat0}, {1.2, 0}],
      Text[Subscript[aIt, "nat0"] // textStyle, {aNat0, 0}, {0, 2.3}]
    },
    (* Regime labels *)
    Graphics @ {
      Text["hot" // textStyle, {Way[aMin, aNatPi], 3/5}],
      Text["warm" // textStyle, {Way[aNatPi, aNat0], 3/4}]
    },
    (* Don't chop labels *)
    ImageMargins -> None,
    PlotRange -> All,
    PlotRangeClipping -> False
  ] // PrettyString["nat" -> "\[Natural]", "pi" -> "\[Pi]"]
] // Ex["bipolar-candidate-asymmetry.pdf"]
