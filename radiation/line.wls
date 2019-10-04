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
(*r_i (inflection radius)*)


(* ::Text:: *)
(*See (r2.32).*)


rInfl = Exp[-2];


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
  With[{a = aHot, phi = \[FormalPhi]},
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
  With[{a = aHot, phi = \[FormalPhi]},
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
  With[{a = aColdHot, phi = \[FormalPhi]},
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
  With[{a = aColdHot, phi = \[FormalPhi]},
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
  With[{a = aCold, phi = \[FormalPhi]},
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
(*Global styles for plots*)


textStyle = Style[#, 18] &;


natStyle = Red;
sharpStyle = Magenta;
flatStyle = Darker[Green];


nonStyle = Directive[Opacity[0.7], LightGray];
termStyle = Pink;
unphysStyle = Black;
inflStyle = Directive[Dashed, Green];


upperStyle = Blue;
lowerStyle = Red;
convexStyle = Black;


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


(* ::Subsection:: *)
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


(* ::Subsection:: *)
(*Series solution for hot regime outer boundaries*)


(* ::Text:: *)
(*See (r2.33) (Page r2-9).*)


With[{r = \[FormalR], xi = \[FormalXi], l = \[FormalCapitalL]},
  Module[{nMax, psiSeries, a, sqrtSeries, phiDerSeries, phiSeries, negLog},
    nMax = 8;
    (* Series for \[Psi](r_\[Sharp] + \[Xi]) *)
    psiSeries = psi[r + xi] + O[xi] ^ nMax;
    (* By construction, \[Psi](r_\[Sharp]) == A *)
    a = psi[r];
    (* Series for sqrt(A^2 - \[Psi]^2(r_\[Sharp] + \[Xi])) *)
    sqrtSeries = Sqrt[a^2 - psiSeries^2];
    (* Series for d\[Phi]/d\[Xi] (see (r2.23a)) *)
    phiDerSeries = psiSeries / ((r + xi) sqrtSeries);
    (* Series for \[Phi] *)
    phiSeries = Integrate[phiDerSeries, xi];
    (* Simplification *)
    phiSeries = phiSeries /. {Log[r] -> -negLog};
    phiSeries =
      Block[{$Assumptions = rNat < r < 1 && 0 < negLog < 4},
        FullSimplify /@ phiSeries
      ];
    phiSeries = phiSeries /. {negLog -> -Log[r]};
    phiSeries = Factor /@ phiSeries;
    phiSeries /. {
      Log[r] -> Subscript[l, "sharp"],
      r -> Subscript[r, "sharp"]
    } // PrettyString["sharp" -> "\[Sharp]"]
  ]
] // Ex["line-hot-traced-outer-series.pdf"]


(* ::Section:: *)
(*Traced boundary plots*)


(* ::Subsection:: *)
(*Hot regime A < A_\[Natural]*)


(* ::Subsubsection:: *)
(*General*)


Module[
 {a,
  phiOffsetValues,
  eps, rMaxAss, rMaxUnphys, rMaxNonAss,
  rMax, rMaxNon, phi
 },
  a = aHot;
  phiOffsetValues = Range[0, 360, 30] Degree // Most;
  eps = 0.1;
  rMaxAss = Association[
    "full" -> 1,
    "zoom_1" -> rInfl,
    "zoom_2" -> rFlat[a]
  ];
  rMaxNonAss = Association[
    "full" -> rSharp[a],
    "zoom_1" -> rInfl (1 + eps),
    "zoom_2" -> rFlat[a] (1 + eps)
  ];
  rMaxUnphys = 1 + eps;
  Table[
    rMax = rMaxAss[zoom];
    rMaxNon = rMaxNonAss[zoom];
    Show[
      EmptyFrame[{-rMax, rMax}, {-rMax, rMax},
        PlotLabel -> BoxedLabel[aIt == a // transStyle[a]]
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
      (* Traced boundaries *)
      Table[
        phi = phiTraHot[id];
        ParametricPlot[
          {
            XYPolar[r, phi[r] + phiOffset],
            XYPolar[r, -phi[r] + phiOffset]
          } // Evaluate, {r, DomainStart[phi], DomainEnd[phi]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {phiOffset, phiOffsetValues}
      , {id, {"outer", "inner"}}],
      (* Radius of inflection *)
      ContourPlot[RPolar[x, y] == rInfl,
        {x, -rInfl, rInfl}, {y, -rInfl, rInfl},
        ContourStyle -> inflStyle
      ]
    ] // Ex["line-hot-traced-" <> zoom <> ".pdf"]
  , {zoom, Keys[rMaxAss]}]
]


(* ::Subsubsection:: *)
(*Convex (generic) without known solution*)


Module[
 {a, rSh,
  nSpike, modNSpike, mod2Pi,
  rSpike, phiSpike, phiCentre, phiHalfWidth,
  rMax, eps, rMaxNon, phi
 },
  a = aHot;
  rSh = rSharp[a];
  nSpike = 3;
  modNSpike = Mod[#, nSpike, 1] &;
  mod2Pi = Mod[#, 2 Pi] &;
  rSpike[1] = rInfl;
  rSpike[2] = Way[rSh, rInfl, 3/4];
  rSpike[3] = Way[rSh, rInfl, 1/4];
  Table[
    phiSpike[n] = phiTraHot["outer"] @ rSpike[n]
  , {n, nSpike}];
  phiCentre[1] = 0;
  phiCentre[2] = 140 Degree;
  phiCentre[3] = 250 Degree;
  phiCentre[n_] := phiCentre[n // modNSpike];
  Table[
    phiHalfWidth[n] =
      phiTraHot["outer"][(1 + 10^-6) rSh] - phiSpike[n]
  , {n, nSpike}];
  phiHalfWidth[n_] := phiHalfWidth[n // modNSpike];
  rMax = rInfl;
  eps = 0.1;
  rMaxNon = rInfl (1 + eps);
  Show[
    EmptyFrame[{-rMax, rMax}, {-rMax, rMax},
      PlotLabel -> BoxedLabel[aIt == a // transStyle[a]]
    ],
    (* Non-viable domain *)
    RegionPlot[vi[a] @ RPolar[x, y] < 0,
      {x, -rMaxNon, rMaxNon}, {y, -rMaxNon, rMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Traced boundaries (spikes) *)
    Table[
      phi[r_] := phiTraHot["outer"][r] - phiSpike[n];
      ParametricPlot[
        {
          XYPolar[r, phiCentre[n] +  phi[r]],
          XYPolar[r, phiCentre[n] -  phi[r]]
        } // Evaluate, {r, rSh, rSpike[n]},
        PlotStyle -> convexStyle
      ]
    , {n, nSpike}],
    (* Traced boundaries (terminal curve portions) *)
    Table[
      ParametricPlot[
        XYPolar[rSh, ph] // Evaluate,
        {
          ph,
          phiCentre[n] + phiHalfWidth[n] // mod2Pi,
          phiCentre[n + 1] - phiHalfWidth[n + 1] // mod2Pi
        },
        PlotStyle -> convexStyle
      ]
    , {n, nSpike}],
    (* Radius of inflection *)
    ContourPlot[RPolar[x, y] == rInfl,
      {x, -rInfl, rInfl}, {y, -rInfl, rInfl},
      ContourStyle -> inflStyle
    ]
  ]
] // Ex["line-hot-traced-convex-generic.pdf"]


(* ::Subsubsection:: *)
(*Convex (generic) with known solution*)


Module[
 {a, rSh,
  nSpike, modNSpike, mod2Pi,
  rSpike, phiSpike, phiCentre, phiHalfWidth,
  rMax, eps, rMaxNon, rBath, phi
 },
  a = aHot;
  rSh = rSharp[a];
  nSpike = 3;
  modNSpike = Mod[#, nSpike, 1] &;
  mod2Pi = Mod[#, 2 Pi] &;
  rSpike[1] = rInfl;
  rSpike[2] = Way[rSh, rInfl, 3/4];
  rSpike[3] = Way[rSh, rInfl, 1/4];
  Table[
    phiSpike[n] = phiTraHot["outer"] @ rSpike[n]
  , {n, nSpike}];
  phiCentre[1] = 0;
  phiCentre[2] = 140 Degree;
  phiCentre[3] = 250 Degree;
  phiCentre[n_] := phiCentre[n // modNSpike];
  Table[
    phiHalfWidth[n] =
      phiTraHot["outer"][(1 + 10^-6) rSh] - phiSpike[n]
  , {n, nSpike}];
  phiHalfWidth[n_] := phiHalfWidth[n // modNSpike];
  rMax = rInfl;
  eps = 0.1;
  rMaxNon = rInfl (1 + eps);
  rBath = rSh / 2;
  Show[
    EmptyFrame[{-rMax, rMax}, {-rMax, rMax},
      PlotLabel -> BoxedLabel[aIt == a // transStyle[a]]
    ],
    (* Non-viable domain *)
    RegionPlot[vi[a] @ RPolar[x, y] < 0,
      {x, -rMaxNon, rMaxNon}, {y, -rMaxNon, rMaxNon},
      BoundaryStyle -> termStyle,
      PlotStyle -> nonStyle
    ],
    (* Known solution *)
    DensityPlot[-Log @ RPolar[x, y],
      {x, -rMax, rMax}, {y, -rMax, rMax},
      ColorFunction -> "TemperatureMap",
      RegionFunction -> Function[{x, y},
        RPolar[x, y] > rBath && Or[
          RPolar[x, y] < rSh,
          Sequence @@ Table[
            phi[r_] := phiTraHot["outer"][r] - phiSpike[n];
            Less[
              phiCentre[n] - phi @ RPolar[x, y],
              PhiPolar[x, y] // If[n == 3, mod2Pi, Identity],
              phiCentre[n] + phi @ RPolar[x, y]
            ]
          , {n, nSpike}]
        ] // Evaluate
      ]
    ],
    (* Traced boundaries (spikes) *)
    Table[
      phi[r_] := phiTraHot["outer"][r] - phiSpike[n];
      ParametricPlot[
        {
          XYPolar[r, phiCentre[n] +  phi[r]],
          XYPolar[r, phiCentre[n] -  phi[r]]
        } // Evaluate, {r, rSh, rSpike[n]},
        PlotStyle -> convexStyle
      ]
    , {n, nSpike}],
    (* Traced boundaries (terminal curve portions) *)
    Table[
      ParametricPlot[
        XYPolar[rSh, ph] // Evaluate,
        {
          ph,
          phiCentre[n] + phiHalfWidth[n] // mod2Pi,
          phiCentre[n + 1] - phiHalfWidth[n + 1] // mod2Pi
        },
        PlotStyle -> convexStyle
      ]
    , {n, nSpike}],
    (* Radius of inflection *)
    ContourPlot[RPolar[x, y] == rInfl,
      {x, -rInfl, rInfl}, {y, -rInfl, rInfl},
      ContourStyle -> inflStyle
    ]
  ]
] // Ex["line-hot-traced-convex-generic-with-known.png"]


(* ::Subsection:: *)
(*Cold-to-hot transition A = A_\[Natural]*)


(* ::Subsubsection:: *)
(*General*)


Module[
 {a,
  phiOffsetValues,
  eps, rMaxAss, rMaxNonAss, rMaxUnphys,
  rMax, rMaxNon, phi
 },
  a = aColdHot;
  phiOffsetValues = Range[0, 360, 30] Degree // Most;
  eps = 0.1;
  rMaxAss = Association[
    "full" -> 1,
    "zoom_1" -> rInfl,
    "zoom_2" -> rFlat[a]
  ];
  rMaxNonAss = Association[
    "full" -> rSharp[a],
    "zoom_1" -> rInfl (1 + eps),
    "zoom_2" -> rFlat[a] (1 + eps)
  ];
  rMaxUnphys = 1 + eps;
  Table[
    rMax = rMaxAss[zoom];
    rMaxNon = rMaxNonAss[zoom];
    Show[
      EmptyFrame[{-rMax, rMax}, {-rMax, rMax},
        PlotLabel -> BoxedLabel[aIt == N[a] // transStyle[a]]
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
      (* Traced boundaries *)
      Table[
        phi = phiTraColdHot[id];
        ParametricPlot[
          {
            XYPolar[r, phi[r] + phiOffset],
            XYPolar[r, -phi[r] + phiOffset]
          } // Evaluate, {r, DomainStart[phi], DomainEnd[phi]},
          PlotStyle -> {upperStyle, lowerStyle},
          RegionFunction -> Function[{x, y},
            Abs[RPolar[x, y] - rNat] > 10^-3 // Evaluate
          ]
        ]
      , {phiOffset, phiOffsetValues}
      , {id, {"outer", "inner"}}],
      (* Radius of inflection *)
      ContourPlot[RPolar[x, y] == rInfl,
        {x, -rInfl, rInfl}, {y, -rInfl, rInfl},
        ContourStyle -> inflStyle
      ],
      (* Transition terminal curve *)
      ContourPlot[RPolar[x, y] == rNat,
        {x, -rMaxNon, rMaxNon}, {y, -rMaxNon, rMaxNon},
        ContourStyle -> termStyle
      ]
    ] // Ex["line-cold_hot-traced-" <> zoom <> ".pdf"]
  , {zoom, Keys[rMaxAss]}]
]


(* ::Subsection:: *)
(*Cold regime A > A_\[Natural]*)


(* ::Subsubsection:: *)
(*General*)


Module[
 {a,
  phiOffsetValues,
  eps, rMaxAss, rMaxUnphys,
  rMax, phi
 },
  a = aCold;
  phiOffsetValues = Range[0, 360, 30] Degree // Most;
  eps = 0.1;
  rMaxAss = Association[
    "full" -> 1,
    "zoom_1" -> rInfl,
    "zoom_2" -> 2 rNat
  ];
  rMaxUnphys = 1 + eps;
  Table[
    rMax = rMaxAss[zoom];
    Show[
      EmptyFrame[{-rMax, rMax}, {-rMax, rMax},
        PlotLabel -> BoxedLabel[aIt == a // transStyle[a]]
      ],
      (* Unphysical domain *)
      RegionPlot[RPolar[x, y] > 1,
        {x, -rMaxUnphys, rMaxUnphys}, {y, -rMaxUnphys, rMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ],
      (* Traced boundaries *)
      Table[
        phi = phiTraCold;
        ParametricPlot[
          {
            XYPolar[r, phi[r] + phiOffset],
            XYPolar[r, -phi[r] + phiOffset]
          } // Evaluate, {r, DomainStart[phi], DomainEnd[phi]},
          PlotStyle -> {upperStyle, lowerStyle}
        ]
      , {phiOffset, phiOffsetValues}],
      (* Radius of inflection *)
      ContourPlot[RPolar[x, y] == rInfl,
        {x, -rInfl, rInfl}, {y, -rInfl, rInfl},
        ContourStyle -> inflStyle
      ]
    ] // Ex["line-cold-traced-" <> zoom <> ".pdf"]
  , {zoom, Keys[rMaxAss]}]
]
