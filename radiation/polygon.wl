(* ::Package:: *)

(* ::Text:: *)
(*See r4 (manuscripts/radiation-4-polygon.pdf).*)


(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< Curvilinear`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "polygon"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Schwarz--Christoffel transformation*)


(* ::Subsubsection:: *)
(*H_n(\[Zeta])*)


(* ::Text:: *)
(*See (r4.8) (Page r4-2).*)


h[n_][zeta_] := Hypergeometric2F1[1/n, 2/n, 1 + 1/n, zeta ^ n];


(* ::Subsubsection:: *)
(*dz/d\[Zeta]*)


(* ::Text:: *)
(*See (r4.11) (Page r4-3).*)


zMapDer[n_][zeta_] := 1 / (h[n][1] (1 - zeta ^ n) ^ (2/n)) // Evaluate;


(* ::Subsubsection:: *)
(*z = z(\[Zeta])*)


(* ::Text:: *)
(*See (r4.12) (Page r4-3).*)


zMap[n_][zeta_] := zeta h[n][zeta] / h[n][1] // Evaluate;


(* ::Subsubsection:: *)
(*d\[Zeta]/dz*)


(* ::Text:: *)
(*See (r4.13) (Page r4-3).*)


zetaMapDer[n_][zeta_] := h[n][1] (1 - zeta ^ n) ^ (2/n) // Evaluate;


(* ::Subsubsection:: *)
(*(Check these)*)


With[{n = \[FormalN], zeta = \[FormalZeta]},
  {
    zMap[n]'[zeta] == zMapDer[n][zeta] // FullSimplify,
    zMapDer[n][zeta] zetaMapDer[n][zeta] == 1
  }
]


(* ::Subsection:: *)
(*Backward Schwarz--Christoffel transformation*)


(* ::Text:: *)
(*Solves z(\[Zeta]) = z by calling FindRoot with initial guess \[Zeta] = z.*)
(*See tests in "Backward transformation" below (zetaTest1 was best).*)


zetaMap[n_][z_] := SeekRoot[zMap[n][#] - z &, z];


(* ::Subsection:: *)
(*Known solution (in \[Zeta]-space)*)


(* ::Subsubsection:: *)
(*G = -log(\[Zeta])*)


(* ::Text:: *)
(*See (r4.14) (Page r4-3).*)


g[zeta_] := -Log[zeta];


(* ::Subsubsection:: *)
(*dG/dz*)


(* ::Text:: *)
(*See (r4.23) (Page r4-4).*)


gDer[n_][zeta_] := -h[n][1] (1 - zeta ^ n) ^ (2/n) / zeta // Evaluate;


(* ::Subsubsection:: *)
(*dG/d\[Zeta]*)


gDerZeta[zeta_] := D[g[zeta], zeta] // Evaluate;


(* ::Subsection:: *)
(*Flux*)


(* ::Subsubsection:: *)
(*F*)


(* ::Text:: *)
(*See (r4.21) (Page r4-4).*)


f[a_][zeta_] := - (Re @ g[zeta]) ^ 4 / a // Evaluate;


(* ::Subsection:: *)
(*Viable domain*)


(* ::Subsubsection:: *)
(*\[CapitalPhi] (viability)*)


(* ::Text:: *)
(*See (r4.22) (Page r4-4).*)


vi[n_][a_][zeta_] :=
  Subtract[
    (Abs @ gDer[n][zeta]) ^ 2,
    (Re @ g[zeta]) ^ 8 / a^2
  ] // Evaluate;


(* ::Subsubsection:: *)
(*\[Psi]*)


(* ::Text:: *)
(*See (r4.24) (Page r4-4).*)


psi[n_][zeta_] := (Re @ g[zeta]) ^ 4 / (Abs @ gDer[n][zeta]);


(* ::Subsubsection:: *)
(*\[Rho]_\[Natural]*)


(* ::Text:: *)
(*The radius \[Rho] = |\[Zeta]| achieving the maximum \[Psi] along arg(\[Zeta]) = \[CurlyPhi].*)
(*See (r4.30) (Page r4-6) and "polygon-psi-algebra.pdf".*)


rhoNat[n_?NumericQ][ph_?NumericQ] :=
  SeekRoot[
    -4 + 8 # ^ n Cos[n ph] + # ^ (2 n) (-4 + Log[#]) - Log[#] &,
    {0, 0.1},
    100
  ]


(* ::Subsubsection:: *)
(*A_\[Natural]*)


(* ::Text:: *)
(*Maximum \[Psi] along arg(\[Zeta]) = \[CurlyPhi].*)


aNat[n_?NumericQ][ph_?NumericQ] := psi[n][rhoNat[n][ph] Exp[I ph]];


(* ::Subsubsection:: *)
(*\[Rho]_\[Sharp]*)


(* ::Text:: *)
(*The larger solution to \[Psi] = A along arg(\[Zeta]) = \[CurlyPhi].*)


rhoSharp[n_?NumericQ][a_?NumericQ][ph_?NumericQ] :=
  Quiet[
    SeekRoot[
      psi[n][# Exp[I ph]] - a &,
      {1, rhoNat[n][ph]},
      100
    ]
  , {Power::infy, Infinity::indet}];


(* ::Subsection:: *)
(*Representative values of A*)


(* ::Text:: *)
(*"Hot" means A < A_\[Natural] along \[CurlyPhi] = 0.*)
(*We take A = A_\[Natural] / 2 as the representative value.*)


aHot[n_] := aNat[n][0] / 2;


(* ::Subsection:: *)
(*Starting points for boundary tracing*)


(* ::Text:: *)
(*We choose points with \[Rho] between \[Rho]_\[Sharp] and 1,*)
(*and \[CurlyPhi] between 0 and 2 \[Pi] / n (since there is n-fold symmetry).*)


(* ::Subsubsection:: *)
(*Hot regime*)


Table[
  startZetaHot[n]["general"] =
    Module[
     {a,
      phMax, phSpacing, phValues,
      rhoMin, rhoMax, rhoBetw, rhoValues
     },
      a = aHot[n];
      (* \[CurlyPhi] values *)
      phMax = 2 Pi / n;
      phSpacing = 30 Degree;
      phValues = UniformRange[0, phMax, phSpacing] // Most;
      (* Build table *)
      Join @@ Table[
        (* \[Rho] values *)
        rhoMin = rhoSharp[n][a][ph];
        rhoMax = 1;
        rhoBetw = Way[rhoMin, rhoMax, 1/5];
        rhoValues = Join[
          Subdivide[rhoMin, rhoBetw, 4] // Most,
          Subdivide[rhoBetw, rhoMax, 4] // Most
        ] // Select[# <= 1 &];
        (* \[Zeta] values *)
        Table[rho Exp[I ph], {rho, rhoValues}]
      , {ph, phValues}]
    ]
, {n, 3, 7}];


(* ::Subsection:: *)
(*Traced boundaries \[Zeta] = \[Zeta](s)*)


(* ::Text:: *)
(*See (r4.38) (Page r4-8).*)


(* ::Subsubsection:: *)
(*Hot regime*)


With[{zeta = \[FormalZeta]},
  Table[
    Module[{a, idList, zetaInitList, viTol},
      (* A *)
      a = aHot[n];
      (* Group names *)
      idList = {"general"};
      (* Solve for traced boundaries *)
      Table[
        zetaInitList = startZetaHot[n][id];
        viTol = vi[n][a] /@ zetaInitList // Min;
        zetaTraHot[id] =
          Table[
            NDSolveValue[
              {
                zeta'[s] == (
                  Divide[
                    I f[a][#] + Sqrt @ vi[n][a][#],
                    gDerZeta[#]
                  ] & @ zeta[s]
                ),
                zeta[0] == zInit,
                WhenEvent[
                  Or[
                    vi[n][a] @ zeta[s] < viTol,
                    Abs @ zeta[s] > 1
                  ],
                  "StopIntegration"
                ]
              }, zeta, {s, -2 Pi, Pi},
              NoExtrapolation
            ]
          , {zInit, zetaInitList}];
      , {id, idList}];
    ]
  , {n, 3, 7}];
];


(* ::Subsection:: *)
(*Geometric regions*)


(* ::Subsubsection:: *)
(*Unit disk (\[Zeta]-space)*)


diskPoints[num_Integer] := diskPoints[num] = (
  SeedRandom[0];
  RandomPoint[Disk[], num]
);


(* ::Subsubsection:: *)
(*Regular n-gon (z-space)*)


poly[n_Integer] := RegularPolygon[{1, 2 Pi / n}, n];


polyVertices[n_Integer] := CirclePoints[{1, 2 Pi / n}, n];


polyPoints[n_Integer, num_Integer] := polyPoints[n, num] = (
  SeedRandom[0];
  RandomPoint[poly[n], num]
)


(* ::Subsection:: *)
(*Italicised symbols*)


aIt = Italicised["A"];
nIt = Italicised["n"];
zIt = Italicised["z"];


(* ::Subsection:: *)
(*Options for exported GIFS*)


gifOpts = Sequence[
  AnimationRepetitions -> Infinity,
  "DisplayDurations" -> 0.5
];


(* ::Subsection:: *)
(*Global styles for plots*)


contStyle = Gray;
nonStyle = Directive[Opacity[0.7], LightGray];
termStyle = Pink;
physStyle = LightGreen;
unphysStyle = Black;


upperStyle = Blue;
lowerStyle = Red;
convexStyle = Black;


pointStyle = PointSize[Large];


(* ::Section:: *)
(*Schwarz--Christoffel transformation*)


(* ::Subsection:: *)
(*Forward transformation: disk (\[Zeta]-space) to polygon (z-space)*)


(* ::Text:: *)
(*See https://mathematica.stackexchange.com/a/159086*)


Table[
  Module[{rSamp = 4, phiSamp = 4},
    Show[
      ParametricPlot[r Exp[I phi] // zMap[n] // ReIm,
        {r, 0, 1}, {phi, 0, 2 Pi},
        Mesh -> {rSamp - 1, n phiSamp - 1},
        PlotPoints -> {rSamp, n phiSamp + 1},
        PlotOptions[Frame] // Evaluate
      ],
      Graphics @ {Directive[Red, pointStyle],
        Point @ polyVertices[n]
      }
    ]
  ] // Ex @ StringJoin["schwarz-christoffel-forward-", ToString[n], ".pdf"]
, {n, 3, 7}]


(* ::Subsection:: *)
(*Backward transformation: polygon (z-space) to disk (\[Zeta]-space)*)


(* ::Subsubsection:: *)
(*Test method 1: call FindRoot (via SeekRoot)*)


zetaTest1[n_][z0_] := SeekRoot[zMap[n][#] - z0 &, z0];


(* ::Subsubsection:: *)
(*Test method 2: integrate d\[Zeta]/dz from 0 to z*)


zetaTest2[n_][z0_] :=
  With[{zeta = \[FormalW]}, (* NOTE: \[FormalZeta] doesn't work *)
    NDSolveValue[
      {
        (* Integrate d\[Zeta]/dz along straight line from z == 0 to z0 *)
        (* but put z == z_0 * t so that t is real *)
        zeta'[t] == zetaMapDer[n] @ zeta[t] * z0,
        zeta[0] == 0
      }, zeta[1], {t, 0, 1}
    ]
  ];


(* ::Subsubsection:: *)
(*Residuals from disk: \[Zeta]_test(z(\[Zeta])) - \[Zeta]*)


Module[{nValues, zetaValues, zetaTestList, resValues},
  nValues = Range[3, 7];
  zetaValues = Complex @@@ diskPoints[100];
  zetaTestList = {zetaTest1, zetaTest2};
  Table[
    resValues =
      Table[
        zetaTest[n] @ zMap[n] @ zeta - zeta
      , {zeta, zetaValues}] // Abs;
    ListPlot[resValues,
      Joined -> True,
      PlotRange -> All
    ]
  , {n, nValues}, {zetaTest, zetaTestList}]
   // TableForm[#, TableHeadings -> {nValues, zetaTestList}] &
] // Ex["schwarz-christoffel-residuals-disk.pdf"]


(* ::Subsubsection:: *)
(*Residuals from polygon: z(\[Zeta]_test(z)) - z*)


Module[{nValues, zetaTestList, zValues, resValues},
  nValues = Range[3, 7];
  zetaTestList = {zetaTest1, zetaTest2};
  Table[
    zValues = Complex @@@ polyPoints[n, 100];
    resValues =
      Table[
        zMap[n] @ zetaTest[n] @ z - z
      , {z, zValues}] // Abs;
    ListPlot[resValues,
      Joined -> True,
      PlotRange -> All
    ]
  , {n, nValues}, {zetaTest, zetaTestList}]
   // TableForm[#, TableHeadings -> {nValues, zetaTestList}] &
] // Ex["schwarz-christoffel-residuals-poly.pdf"]


(* ::Subsubsection:: *)
(*Conclusion*)


(* ::Text:: *)
(*zetaTest1 has much smaller residuals.*)


(* ::Subsubsection:: *)
(*Timings*)


Module[{nValues, zetaTestList, zValues, timingList},
  nValues = Range[3, 7];
  zetaTestList = {zetaTest1, zetaTest2};
  timingList =
    Transpose @ Table[
      zValues = Complex @@@ polyPoints[n, 1000];
      Table[
        zMap[n] @ zetaTest[n] @ z - z
      , {z, zValues}] // AbsoluteTiming // First
    , {n, nValues}, {zetaTest, zetaTestList}];
  ListPlot[timingList,
    AxesLabel -> {Italicised["n"], "Timing"},
    Joined -> True,
    PlotLegends -> zetaTestList,
    PlotRange -> All,
    PlotOptions[Axes] // Evaluate
  ]
] // Ex["schwarz-christoffel-timings.pdf"]


(* ::Subsubsection:: *)
(*Conclusion*)


(* ::Text:: *)
(*zetaTest1 is faster (and the difference increases with n).*)


(* ::Section:: *)
(*\[Psi]*)


(* ::Subsection:: *)
(*Algebra*)


With[{n = \[FormalN], rho = \[FormalRho], ph = \[FormalCurlyPhi]},
  Block[{$Assumptions = n > 2 && 0 < rho < 1 && -Pi < ph < Pi},
    Module[{psiExpr, psiDerExpr},
      psiExpr = psi[n][rho Exp[I ph]] // ComplexExpand // FullSimplify;
      psiDerExpr = psiExpr // D[#, rho] & // FullSimplify;
      {
        {"psi", "(r4.29)", psiExpr},
        {"dpsi/dr", "(r4.30)", psiDerExpr}
      } /. {Gamma -> "\[CapitalGamma]"}
        // PrettyString["d" -> "\[PartialD]", "psi" -> "\[Psi]"]
        // TableForm
    ]
  ]
] // Ex["polygon-psi-algebra.pdf"]


(* ::Subsection:: *)
(*\[Psi] plot*)


Module[{phValues},
  Table[
    (* \[CurlyPhi] values *)
    phValues = Subdivide[0, Pi / n, 4];
    (* Plot *)
    Plot[
      Table[
        psi[n][rho Exp[I ph]]
      , {ph, phValues}] // Evaluate,
      {rho, 0, 1},
      AxesLabel -> {"rho", "psi"},
      ImageSize -> 360,
      PlotLabel -> nIt == n,
      PlotLegends -> LineLegend[
        phValues,
        LegendLabel -> "ph"
      ],
      PlotRange -> Full,
      PlotOptions[Axes] // Evaluate
    ] // PrettyString[
      "rho" -> "\[Rho]",
      "ph" -> "\[CurlyPhi]",
      "psi" -> "\[Psi]"
    ] // Ex @ StringJoin["polygon-psi-", ToString[n], ".pdf"]
  , {n, 3, 7}]
]


(* ::Subsection:: *)
(*\[Rho]_\[Natural] angular dependence plot*)


Table[
  Plot[rhoNat[n][ph], {ph, 0, 2 Pi},
    AxesLabel -> {"ph", Subscript["rho", "Nat"]},
    ImageSize -> 360,
    PlotLabel -> Column[
      {
        nIt == n,
        Equal[
          Row[{Max, Min}, "/"] - 1,
          rhoNat[n][0] / rhoNat[n][Pi / n] - 1
        ]
      },
      Alignment -> Center
    ],
    PlotRange -> Full,
    PlotOptions[Axes] // Evaluate
  ] // PrettyString[
    "rho" -> "\[Rho]",
    "ph" -> "\[CurlyPhi]",
    "Nat" -> "\[Natural]"
  ] // Ex @ StringJoin["polygon-rho-nat-", ToString[n], ".pdf"]
, {n, 3, 7}]


(* ::Subsection:: *)
(*A_\[Natural] angular dependence plot*)


Table[
  Plot[aNat[n][ph], {ph, 0, 2 Pi},
    AxesLabel -> {"ph", Subscript[aIt, "Nat"]},
    ImageSize -> 360,
    PlotLabel -> Column[
      {
        nIt == n,
        Equal[
          Row[{Max, Min}, "/"] - 1,
          aNat[n][0] / aNat[n][Pi / n] - 1
        ]
      },
      Alignment -> Center
    ],
    PlotRange -> Full,
    PlotOptions[Axes] // Evaluate
  ] // PrettyString[
    "ph" -> "\[CurlyPhi]",
    "Nat" -> "\[Natural]"
  ] // Ex @ StringJoin["polygon-a-nat-", ToString[n], ".pdf"]
, {n, 3, 7}]


(* ::Section:: *)
(*Viable domain*)


(* ::Subsection:: *)
(*Animations (\[Zeta]-space)*)


Table[
  Module[
   {ph,
    aN, aMin, aMax, aStep, aValues,
    eps, rhoMax, rhoMaxUnphys, rhoMaxNon
   },
   (* Azimuthal angle in \[Zeta]-space *)
   ph = 0;
   (* Values of A *)
   aN = aNat[n][ph];
   aMin = 1/10;
   aMax = Ceiling[aN];
   aStep = 1/10;
   aValues = Range[aMin, aMax, aStep];
   aValues = Append[aValues, aN] // Sort[#, Less] &;
   (* Animation *)
   Table[
     eps = 0.1;
     rhoMax = 1;
     rhoMaxUnphys = 1 + eps;
     rhoMaxNon = If[a < aN, rhoSharp[n][a][ph], rhoNat[n][ph]];
     Show[
       EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
         FrameLabel -> {Re["zeta"], Im["zeta"]},
         ImageSize -> 360,
         PlotLabel -> BoxedLabel[aIt == N[a]]
       ] // PrettyString["zeta" -> "\[Zeta]"],
       (* Unphysical domain *)
       RegionPlot[RPolar[reZeta, imZeta] > 1,
         {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
         {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
         BoundaryStyle -> None,
         PlotStyle -> unphysStyle
       ],
       (* Non-viable domain *)
       RegionPlot[vi[n][a][reZeta + I imZeta] < 0,
         {reZeta, -rhoMaxNon, rhoMaxNon},
         {imZeta, -rhoMaxNon, rhoMaxNon},
         BoundaryStyle -> termStyle,
         PlotStyle -> nonStyle
       ]
     ]
   , {a, aValues}]
  ] // Ex[
    StringJoin["polygon_zeta-viable-full-", ToString[n],".gif"],
    gifOpts
  ]
, {n, 3, 7}]


Table[
  Module[
   {ph,
    aN, aStep, aMin, aMax, aValues,
    eps, rhoN, rhoMax, rhoMaxUnphys, rhoMaxNon
   },
   (* Azimuthal angle in \[Zeta]-space *)
   ph = 0;
   (* Values of A *)
   aN = aNat[n][ph];
   aStep = 1/20;
   aMin = Floor[aN - 5 aStep, aStep];
   aMax = Ceiling[aN, aStep];
   aValues = Range[aMin, aMax, aStep];
   aValues = Append[aValues, aN] // Sort[#, Less] &;
   (* Animation *)
   Table[
     eps = 0.1;
     rhoN = rhoNat[n][ph];
     rhoMax = 2 rhoN;
     rhoMaxUnphys = 1 + eps;
     rhoMaxNon = If[a < aN, rhoSharp[n][a][ph], rhoN];
     Show[
       EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
         FrameLabel -> {Re["zeta"], Im["zeta"]},
         ImageSize -> 360,
         PlotLabel -> BoxedLabel[aIt == N[a]]
       ] // PrettyString["zeta" -> "\[Zeta]"],
       (* Non-viable domain *)
       RegionPlot[vi[n][a][reZeta + I imZeta] < 0,
         {reZeta, -rhoMaxNon, rhoMaxNon},
         {imZeta, -rhoMaxNon, rhoMaxNon},
         BoundaryStyle -> termStyle,
         PlotStyle -> nonStyle
       ]
     ]
   , {a, aValues}]
  ] // Ex[
    StringJoin["polygon_zeta-viable-zoom-", ToString[n],".gif"],
    gifOpts
  ]
, {n, 3, 7}]


(* ::Section:: *)
(*Traced boundary plots*)


(* ::Subsection:: *)
(*Hot regime*)


(* ::Subsubsection:: *)
(*Starting points (\[Zeta]-space)*)


Table[
  Module[
   {rSamp = 4, phiSamp = 4,
    a, idList,
    eps, rhoMax, rhoMaxUnphys, rhoMaxNon
   },
    a = aHot[n];
    (* Group names *)
    idList = {"general"};
    (* Plot ranges *)
    eps = 0.1;
    rhoMax = 1;
    rhoMaxUnphys = 1 + eps;
    rhoMaxNon = rhoSharp[n][a][0];
    (* Plot *)
    Show[
      EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
        FrameLabel -> {Re["zeta"], Im["zeta"]},
        ImageSize -> 360,
        PlotLabel -> BoxedLabel @ Row[
          {aIt == N[aHot[n]], "hot" // ""},
          "  "
        ]
      ] // PrettyString["zeta" -> "\[Zeta]"],
      (* Unphysical domain *)
      RegionPlot[RPolar[reZeta, imZeta] > 1,
        {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
        {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ],
      (* Non-viable domain *)
      RegionPlot[vi[n][a][reZeta + I imZeta] < 0,
        {reZeta, -rhoMaxNon, rhoMaxNon},
        {imZeta, -rhoMaxNon, rhoMaxNon},
        BoundaryStyle -> termStyle,
        PlotStyle -> nonStyle
      ],
      (* Starting points *)
      ListPlot[
        Table[startZetaHot[n][id] // ReIm, {id, idList}],
        LabelingFunction -> Function @ Placed[
          #2[[2]],
          Center
        ],
        PlotLegends -> idList,
        PlotStyle -> Directive[Opacity[0.7], pointStyle]
      ]
    ]
  ] // Ex @ StringJoin["polygon_zeta-hot-traced-starting-", ToString[n],".pdf"]
, {n, 3, 7}]


(* ::Subsubsection:: *)
(*Starting points (z-space)*)


Table[
  Module[{rSamp = 4, phiSamp = 4, idList},
    (* Group names *)
    idList = {"general"};
    (* Plot *)
    Show[
      (* Domain with equipotentials (contours) *)
      ParametricPlot[r Exp[I phi] // zMap[n] // ReIm,
        {r, 0, 1}, {phi, 0, 2 Pi},
        BoundaryStyle -> unphysStyle,
        ImageSize -> 480,
        Mesh -> {rSamp - 1, n phiSamp - 1},
        MeshStyle -> contStyle,
        PlotLabel -> BoxedLabel @ Row[
          {aIt == N[aHot[n]], "hot" // ""},
          "  "
        ],
        PlotPoints -> {rSamp, n phiSamp + 1},
        PlotStyle -> physStyle,
        PlotOptions[Frame] // Evaluate
      ],
      (* Starting points *)
      ListPlot[
        Table[startZetaHot[n][id] // zMap[n] // ReIm, {id, idList}],
        LabelingFunction -> Function @ Placed[
          #2[[2]],
          Center
        ],
        PlotLegends -> idList,
        PlotStyle -> Directive[Opacity[0.7], pointStyle]
      ]
    ]
  ] // Ex @ StringJoin["polygon_z-hot-traced-starting-", ToString[n],".pdf"]
, {n, 3, 7}]


(* ::Subsubsection:: *)
(*General (\[Zeta]-space)*)


Table[
  Module[
   {rSamp = 4, phiSamp = 4,
    a, idList,
    eps, rhoMax, rhoMaxUnphys, rhoMaxNon
   },
    a = aHot[n];
    (* Group names *)
    idList = {"general"};
    (* Plot ranges *)
    eps = 0.1;
    rhoMax = 1;
    rhoMaxUnphys = 1 + eps;
    rhoMaxNon = rhoSharp[n][a][0];
    (* Plot *)
    Show[
      EmptyFrame[{-rhoMax, rhoMax}, {-rhoMax, rhoMax},
        FrameLabel -> {Re["zeta"], Im["zeta"]},
        ImageSize -> 360,
        PlotLabel -> BoxedLabel @ Row[
          {aIt == N[aHot[n]], "hot" // ""},
          "  "
        ]
      ] // PrettyString["zeta" -> "\[Zeta]"],
      (* Unphysical domain *)
      RegionPlot[RPolar[reZeta, imZeta] > 1,
        {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
        {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ],
      (* Non-viable domain *)
      RegionPlot[vi[n][a][reZeta + I imZeta] < 0,
        {reZeta, -rhoMaxNon, rhoMaxNon},
        {imZeta, -rhoMaxNon, rhoMaxNon},
        BoundaryStyle -> termStyle,
        PlotStyle -> nonStyle
      ],
      (* Traced boundaries *)
      Table[
        Table[
          Table[
            ParametricPlot[
              zeta[s] Exp[I 2 Pi k / n]
                // {#, Conjugate[#]} &
                // ReIm
                // Evaluate,
              {s, DomainStart[zeta], DomainEnd[zeta]},
              PlotStyle -> {upperStyle, lowerStyle}
            ]
          , {k, 0, n - 1}]
        , {zeta, zetaTraHot[id]}]
      , {id, idList}]
    ]
  ] // Ex @ StringJoin["polygon_zeta-hot-traced-full-", ToString[n],".pdf"]
, {n, 3, 7}]


(* ::Subsubsection:: *)
(*General (z-space)*)


Table[
  Module[{rSamp = 4, phiSamp = 4, idList},
    (* Group names *)
    idList = {"general"};
    (* Plot *)
    Show[
      (* Domain with equipotentials (contours) *)
      ParametricPlot[r Exp[I phi] // zMap[n] // ReIm,
        {r, 0, 1}, {phi, 0, 2 Pi},
        BoundaryStyle -> unphysStyle,
        ImageSize -> 480,
        Mesh -> {rSamp - 1, n phiSamp - 1},
        MeshStyle -> contStyle,
        PlotLabel -> BoxedLabel @ Row[
          {aIt == N[aHot[n]], "hot" // ""},
          "  "
        ],
        PlotPoints -> {rSamp, n phiSamp + 1},
        PlotStyle -> physStyle,
        PlotOptions[Frame] // Evaluate
      ],
      (* Traced boundaries *)
      Table[
        Table[
          Table[
            ParametricPlot[
              zeta[s] Exp[I 2 Pi k / n]
                // zMap[n]
                // {#, Conjugate[#]} &
                // ReIm
                // Evaluate,
              {s, DomainStart[zeta], DomainEnd[zeta]},
              PlotStyle -> {upperStyle, lowerStyle}
            ]
          , {k, 0, n - 1}]
        , {zeta, zetaTraHot[id]}]
      , {id, idList}]
    ]
  ] // Ex @ StringJoin["polygon_z-hot-traced-full-", ToString[n],".pdf"]
, {n, 3, 7}]
