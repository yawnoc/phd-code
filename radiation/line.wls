#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*See r2 (manuscripts/radiation-2-line.pdf).*)


(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< Curvilinear`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "line"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*\[Psi] = \[Psi](r)*)


(* ::Text:: *)
(*See (r2.16).*)


psi[r_] := r Log[r]^4;


(* ::Subsection:: *)
(*r_\[Natural] (transition terminal radius)*)


(* ::Text:: *)
(*See (r2.19).*)


rNat = Exp[-4]


(* ::Subsection:: *)
(*A_\[Natural] (transition dimensionless group)*)


(* ::Text:: *)
(*See (r2.20).*)


aNat = psi[rNat]


(* ::Subsection:: *)
(*r_\[Flat] (inner terminal radius)*)


(* InverseFunction can't invert psi ... *)(*
InverseFunction @ Function[{r},
  ConditionalExpression[
    psi[r],
    0 < r < rNat
  ] // Evaluate
]*)


(* ... so spot pattern from special values in Solve: *)(*
With[{r = \[FormalR]},
  Table[
    Select[
      r /. Solve[psi[r] == a, r, Reals],
      0 < # < rNat &
    ]
  , {a, {1, E, 3, Pi}}]
]*)


rFlat[a_] := Exp[4 ProductLog[-1, -a^(1/4) / 4]];


(* Check this: *)(*
With[{a = \[FormalCapitalA]},
  Block[{$Assumptions = 0 < a < aNat},
    psi @ rFlat[a] == a // FullSimplify
  ]
]*)


(* ::Subsection:: *)
(*r_\[Sharp] (outer terminal radius)*)


(* InverseFunction can't invert psi ... *)(*
InverseFunction @ Function[{r},
  ConditionalExpression[
    psi[r],
    rNat < r < 1
  ] // Evaluate
]*)


(* ... so spot pattern from special values in Solve: *)(*
With[{r = \[FormalR]},
  Table[
    Select[
      r /. Solve[psi[r] == a, r, Reals],
      rNat < # < 1 &
    ]
  , {a, {1, E, 3, Pi}}]
]*)


rSharp[a_] := Exp[4 ProductLog[-a^(1/4) / 4]];


(* Check this: *)(*
With[{a = \[FormalCapitalA]},
  Block[{$Assumptions = 0 < a < aNat},
    psi @ rSharp[a] == a // FullSimplify
  ]
]*)


(* ::Subsection:: *)
(*\[CapitalPhi] (viability)*)


(* phi I am reserving for azimuthal angle *)
vi[a_] := Function[{r},
  (a^2 - psi[r]^2) / (a^2 r^2) // Evaluate
] // Evaluate;


(* ::Subsection:: *)
(*Representative values of A*)


aHot = 3;
aColdHot = aNat;
aCold = 5;
aHot < aColdHot < aCold


(* ::Subsection:: *)
(*Traced boundary \[Phi] = \[Phi](r)*)


phiTraDer[a_] := Function[{r},
  -psi[r] / (r Sqrt[a^2 - psi[r]^2]) // Evaluate
] // Evaluate;


(* ::Subsubsection:: *)
(*Hot regime A < A_\[Natural]*)


phiTraHot["outer"] =
  With[
   {a = aHot,
    eps = 10^-6,
    phi = \[FormalPhi]
   },
    Module[{rMin, rMax, rInit},
      rMin = rSharp[a];
      rMax = 1;
      rInit = Way[rMin, rMax];
      NDSolveValue[
        {phi'[r] == phiTraDer[a][r], phi[rInit] == 0},
        phi, {r, rMin, rMax},
        NoExtrapolation
      ]
    ]
  ];


phiTraHot["inner"] =
  With[
   {a = aHot,
    eps = 10^-6,
    phi = \[FormalPhi]
   },
    Module[{rMin, rMax, rInit},
      rMin = 0;
      rMax = rFlat[a];
      rInit = Way[rMin, rMax];
      NDSolveValue[
        {phi'[r] == phiTraDer[a][r], phi[rInit] == 0},
        phi, {r, rMin, rMax},
        NoExtrapolation
      ]
    ]
  ];


(* ::Subsubsection:: *)
(*Cold-to-hot transition A = A_\[Natural]*)


phiTraColdHot["outer"] =
  With[
   {a = aColdHot,
    eps = 10^-6,
    phi = \[FormalPhi]
   },
    Module[{rMin, rMax, rInit},
      rMin = rNat;
      rMax = 1;
      rInit = Way[rMin, rMax];
      NDSolveValue[
        {phi'[r] == phiTraDer[a][r], phi[rInit] == 0},
        phi, {r, rMin, rMax},
        NoExtrapolation
      ]
    ]
  ];


phiTraColdHot["inner"] =
  With[
   {a = aColdHot,
    eps = 10^-6,
    phi = \[FormalPhi]
   },
    Module[{rMin, rMax, rInit},
      rMin = 0;
      rMax = rNat;
      rInit = Way[rMin, rMax];
      NDSolveValue[
        {phi'[r] == phiTraDer[a][r], phi[rInit] == 0},
        phi, {r, rMin, rMax},
        NoExtrapolation
      ]
    ]
  ];


(* ::Subsubsection:: *)
(*Cold regime A > A_\[Natural]*)


phiTraCold =
  With[
   {a = aCold,
    eps = 10^-6,
    phi = \[FormalPhi]
   },
    Module[{rMin, rMax, rInit},
      rMin = 0;
      rMax = 1;
      rInit = Way[rMin, rMax];
      NDSolveValue[
        {phi'[r] == phiTraDer[a][r], phi[rInit] == 0},
        phi, {r, rMin, rMax},
        NoExtrapolation
      ]
    ]
  ];


(* ::Subsection:: *)
(*Italicised symbols*)


aIt = Italicised["A"];
rIt = Italicised["r"];


(* ::Subsection:: *)
(*Options for exported GIFS*)


gifOpts = Sequence[
  AnimationRepetitions -> Infinity,
  "DisplayDurations" -> 0.5
];


(* ::Subsection:: *)
(*Global styles*)


textStyle = Style[#, 18] &;


natStyle = Red;
sharpStyle = Magenta;
flatStyle = Darker[Green];


nonStyle = Directive[Opacity[0.7], LightGray];
termStyle = Pink;
unphysStyle = Black;


upperStyle = Blue;
lowerStyle = Red;


transStyle = Function[{a}, If[a == aNat, Style[#, Red] &, Identity]];


(* ::Section:: *)
(*\[Psi]*)


(* ::Subsection:: *)
(*Properties*)


With[{r = \[FormalR]},
  {
    {"psi", psi[r]},
    {"psi(0)", Limit[psi[r], r -> 0]},
    {"psi(1)", psi[1]},
    {"psi'", psi'[r] // FullSimplify},
    {"psi'(0)", Limit[psi'[r], r -> 0]},
    {"psi'(rNat)", psi'[rNat]},
    {"psi'(1)", psi'[1]},
    {"psi''", psi''[r] // FullSimplify},
    {"psi''(rNat)", psi''[rNat]},
    {"psi''(rNat)", psi''[rNat] // N},
    {"rNat", rNat},
    {"rNat", rNat // N},
    {"ANat", aNat},
    {"ANat", aNat // N}
  } // PrettyString["Nat" -> "_\[Natural]", "psi" -> "\[Psi]"]
    // TableForm
] // Ex["line-psi-properties.pdf"]


(* ::Subsection:: *)
(*Plot*)


Show[
  (* psi *)
  Plot[psi[r], {r, 0, 1},
    AxesLabel -> {rIt, aIt},
    PlotRange -> All,
    PlotRangeClipping -> False,
    PlotStyle -> Blue,
    PlotOptions[Axes] // Evaluate
  ],
  (* rNat *)
  Graphics @ {Directive[Dashed, natStyle],
    Line @ {{rNat, 0}, {rNat, aNat}}
  },
  Graphics @ {natStyle,
    Text[
      Subscript[rIt, "Nat"] // textStyle,
      {rNat, 0},
      {-0.7, 0.8}
    ]
  },
  (* aNat *)
  Graphics @ {Directive[Dashed, natStyle],
    Line @ {{rNat, aNat}, {0, aNat}}
  },
  Graphics @ {natStyle,
    Text[
      Subscript[aIt, "Nat"] // textStyle,
      {0, aNat},
      {1.3, 0}
    ]
  }
] // PrettyString["Nat" -> "\[Natural]"] // Ex["line-psi.pdf"]


(* ::Section:: *)
(*Viable domain*)


(* ::Subsection:: *)
(*Plot (radial space)*)


Module[{eps},
  eps = 0.07;
  Show[
    EmptyFrame[{0, (1 + eps) aNat}, {0, 1 + eps},
      AspectRatio -> 1,
      FrameLabel -> {aIt, rIt}
    ],
    (* Viable domain *)
    Graphics @ {
      Text[
        "Viable domain" // textStyle,
        {(1 + eps) aNat / 2, (1 + eps) / 2}
      ]
    },
    (* Terminal radii *)
    Plot[{rSharp[a], rFlat[a]}, {a, 0, aNat},
    (* Non-viable domain *)
      Filling -> {1 -> {2}},
      FillingStyle -> nonStyle,
      PlotLegends -> textStyle /@ {
        rIt == Subscript[rIt, "Sharp"],
        rIt == Subscript[rIt, "Flat"]
      },
      PlotRange -> Full,
      PlotStyle -> {sharpStyle, flatStyle}
    ],
    Graphics @ {
      Text[
        "Non-viable domain" // textStyle,
        {(1 + eps) aNat / 3, (1 + eps) / 20}
      ]
    },
    (* Unphysical domain *)
    Plot[1, {a, 0, (1 + 2 eps) aNat},
      Filling -> Top,
      FillingStyle -> unphysStyle,
      PlotStyle -> unphysStyle
    ],
    Graphics @ {White,
      Text[
        "Unphysical domain" // textStyle,
        {(1 + eps) aNat / 2, 1 + eps / 2}
      ]
    }
  ] // PrettyString["Sharp" -> "\[Sharp]", "Flat" -> "\[Flat]"]
] // Ex["line-viable-r.pdf"]


(* ::Subsection:: *)
(*Animations*)


Module[
 {aMin, aMax, aStep,
  aValues,
  eps, rMax, rMaxUnphys, rMaxNon
 },
  (* Values of A *)
  aMin = 1/5;
  aMax = Ceiling[aNat];
  aStep = 1/5;
  aValues = Range[aMin, aMax, aStep];
  aValues = Append[aValues, aNat] // Sort[#, Less] &;
  (* Animation *)
  Table[
    eps = 0.1;
    rMax = 1;
    rMaxUnphys = 1 + eps;
    rMaxNon = If[a < aNat, rSharp[a], rNat];
    Show[
      EmptyFrame[{-rMax, rMax}, {-rMax, rMax},
        PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
      ],
      (* Unphysical domain *)
      RegionPlot[RPolar[x, y] > 1,
        {x, -rMaxUnphys, rMaxUnphys}, {y, -rMaxUnphys, rMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ],
      (* Non-viable domain *)
      RegionPlot[vi[a] @ RPolar[x, y] < 0,
        {x, -rMaxNon, rMaxNon}, {y, -rMaxNon, rMaxNon},
        BoundaryStyle -> termStyle,
        PlotStyle -> nonStyle
      ],
      (* Transition terminal curve *)
      If[a == aNat,
        ContourPlot[RPolar[x, y] == rNat,
          {x, -rMaxNon, rMaxNon}, {y, -rMaxNon, rMaxNon},
          ContourStyle -> termStyle
        ],
        {}
      ]
    ]
  , {a, aValues}]
] // Ex["line-viable-full.gif", gifOpts]


Module[
 {aMin, aMax, aStep,
  aValues,
  eps, rMax, rMaxUnphys, rMaxNon
 },
  (* Values of A *)
  aMin = Floor[aNat];
  aMax = Ceiling[aNat];
  aStep = 1/10;
  aValues = Range[aMin, aMax, aStep];
  aValues = Append[aValues, aNat] // Sort[#, Less] &;
  (* Animation *)
  Table[
    eps = 0.1;
    rMax = 2 rNat;
    rMaxUnphys = 1 + eps;
    rMaxNon = If[a < aNat, rSharp[a], rNat];
    Show[
      EmptyFrame[{-rMax, rMax}, {-rMax, rMax},
        PlotLabel -> BoxedLabel @ transStyle[a][aIt == N[a]]
      ],
      (* Unphysical domain *)
      RegionPlot[RPolar[x, y] > 1,
        {x, -rMaxUnphys, rMaxUnphys}, {y, -rMaxUnphys, rMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ],
      (* Non-viable domain *)
      RegionPlot[vi[a] @ RPolar[x, y] < 0,
        {x, -rMaxNon, rMaxNon}, {y, -rMaxNon, rMaxNon},
        BoundaryStyle -> termStyle,
        PlotStyle -> nonStyle
      ],
      (* Transition terminal curve *)
      If[a == aNat,
        ContourPlot[RPolar[x, y] == rNat,
          {x, -rMaxNon, rMaxNon}, {y, -rMaxNon, rMaxNon},
          ContourStyle -> termStyle
        ],
        {}
      ]
    ]
  , {a, aValues}]
] // Ex["line-viable-zoom.gif", gifOpts]


(* ::Section:: *)
(*Traced boundary algebra*)


(* ::Subsubsection:: *)
(*Curvature and radius of inflection*)


(* ::Text:: *)
(*See Page r2-8.*)


With[{r = \[FormalR], phi = \[FormalPhi], a = \[FormalCapitalA]},
  Module[
   {ar, aphi, d, dSimp, vel, acc, comp, cur,
    phiDer, phiDer2, curEval, rInfl,
    brevity
   },
    (* Local orthogonal basis *)
    {ar, aphi};
    (* Abbreviation for r-derivative *)
    d = Dt[#, r, Constants -> a] &;
    (* Simplification rules *)
    dSimp = {
      d[ar] -> phi'[r] aphi,
      d[aphi] -> -phi'[r] ar
    };
    (* Velocity and acceleration *)
    vel = ar + r phi'[r] aphi;
    acc = d[vel] /. dSimp;
    (* Components for local (r, \[Phi], z) basis *)
    comp = {
      ar -> {1, 0, 0},
      aphi -> {0, 1, 0}
    };
    {vel, acc} = {vel, acc} /. comp;
    (* Curvature *)
    cur = Last @ Cross[vel, acc];
    (* \[Phi]' and \[Phi]'' *)
    phiDer = phiTraDer[a][r];
    phiDer2 = d[phiDer] // FullSimplify;
    (* Curvature evaluated *)
    curEval = cur /. {
      phi'[r] -> phiDer,
      phi''[r] -> phiDer2
    } // FullSimplify;
    (* Radius of inflection *)
    rInfl = r /. Solve[curEval == 0 && 0 < r < 1, r];
    (* For brevity *)
    brevity = {Log[r] -> \[FormalCapitalL]};
    (* Table of results *)
    {
      {"Curvature", cur},
      {"phi'", phiDer},
      {"phi''", phiDer2},
      {"Curvature", curEval},
      {"r_i", rInfl}
    } /. brevity
      // PrettyString["phi" -> "\[Phi]"]
      // TableForm
  ]
] // Ex["line-curvature-inflection.pdf"]


(* ::Section:: *)
(*Traced boundary plots*)


(* ::Subsubsection:: *)
(*Hot regime A < A_\[Natural]*)
