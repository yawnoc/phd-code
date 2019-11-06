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
      rhoMin, rhoMax, rho,
      phMax, phSpacing, phValues
     },
      a = aHot[n];
      (* \[Rho] *)
      rhoMin = rhoSharp[n][a][0];
      rhoMax = 1;
      rho = Way[rhoMin, rhoMax, 1/2];
      (* \[CurlyPhi] values *)
      phMax = 2 Pi / n;
      phSpacing = 30 Degree;
      phValues = UniformRange[0, phMax, phSpacing] // Most;
      (* Build list of values *)
      Table[rho Exp[I ph], {ph, phValues}]
    ]
, {n, 3, 5}];


Table[
  startZetaHot[n]["hyperbolic"] =
    Module[{a, rhoSh},
      a = aHot[n];
      (* Build list of values *)
      (* (only need singleton due to rotational symmetry) *)
      {rhoSharp[n][a][0]}
    ]
, {n, 3, 5}];


(* ::Subsection:: *)
(*Traced boundaries \[Zeta] = \[Zeta](s)*)


(* ::Subsubsection:: *)
(*\[Zeta]' = d\[Zeta]/ds*)


(* ::Text:: *)
(*See (r4.38) (Page r4-8).*)


zetaVel[n_][a_][zeta_] :=
  Divide[
    I f[a][zeta] + Sqrt @ vi[n][a][zeta],
    gDerZeta[zeta]
  ] // Evaluate;


(* ::Subsubsection:: *)
(*Solver for traced boundaries*)


(* ::Text:: *)
(*See (r4.38) (Page r4-8).*)


With[{zeta = \[FormalZeta]},
  zetaTrace[n_][a_][zetaInit_] :=
    NDSolveValue[
      {
        zeta'[s] == zetaVel[n][a] @ zeta[s],
        zeta[0] == zetaInit,
        WhenEvent[
          vi[n][a] @ zeta[s] < 0 || Abs @ zeta[s] > 1,
          "StopIntegration"
        ]
      }, zeta, {s, -Pi, Pi},
      NoExtrapolation
    ]
];


(* ::Subsubsection:: *)
(*Solver for hyperbolic traced boundaries*)


zetaTraceHyperbolic[n_][a_] :=
  Module[{rhoClearance},
    rhoClearance = 10^-8;
    zetaTrace[n][a] @ (rhoSharp[n][a][0] + rhoClearance)
  ];


(* ::Subsubsection:: *)
(*Hot regime*)


With[{zeta = \[FormalZeta]},
  Table[
    Module[{a, idList, zetaInitList},
      (* A *)
      a = aHot[n];
      (* Group names *)
      idList = {"general", "hyperbolic"};
      (* Solve for traced boundaries *)
      Table[
        zetaInitList = startZetaHot[n][id];
        zetaTraHot[n][id] =
          Table[
            zetaTrace[n][a][zetaInit]
          , {zetaInit, zetaInitList}];
      , {id, idList}];
    ]
  , {n, 3, 5}];
];


(* ::Subsection:: *)
(*Traced boundary curvature*)


(* ::Subsubsection:: *)
(*\[CurlyPhi](\[Rho] = 1)*)


phEnd[n_][a_] :=
  Module[{zeta},
    zeta = zetaTraceHyperbolic[n][a];
    zeta @ DomainStart[zeta] // Conjugate // Arg
  ];


(* ::Subsubsection:: *)
(*A_m (meeting dimensionless group)*)


(* ::Text:: *)
(*This is the smallest A > 0 for which \[CurlyPhi] = \[Pi]/n at \[Rho] = 1.*)
(*See Page r4-10.*)


(* Compute A_m using the bisection algorithm *)
aMeet = Module[
 {dest, nValues,
  aMin, aMax,
  a, num,
  aAss, numAss
 },
  (* (This is not slow, nevertheless compute once and store.) *)
  (* (Delete the file manually to compute from scratch.) *)
  dest = "polygon-a-meet.txt";
  nValues = Range[3, 5];
  If[FileExistsQ[dest],
    (* If already stored, import *)
    {aAss, numAss} = Import[dest] // Uncompress,
    (* Otherwise compute and export for next time *)
    aAss = Association[];
    numAss = Association[];
    aMin = 5/10;
    aMax = 12/10;
    (* NOTE: these may not work for n > 5 *)
    Table[
      {a, num} = SeekRootBisection[
        phEnd[n][#] - Pi / n &,
        {aMin, aMax},
        "ReturnIterations" -> True
      ];
      aAss = aAss // Append[n -> a];
      numAss = numAss // Append[n -> num];
    , {n, nValues}];
    {aAss, numAss} // Compress // Ex[dest]
  ];
  (* Print iterations used *)
  Table[
    Print[
      "n == {n} bisection algorithm: {num} iterations"
        // PrettyString[
          "{n}" -> ToString[n],
          "{num}" -> ToString @ numAss[n]
        ]
    ];
  , {n, nValues}];
  (* Return A_m values *)
  aAss
];
aMeet // N


(* ::Subsubsection:: *)
(*Solver for candidate traced boundaries*)


(* ::Text:: *)
(*This is essentially the same as zetaTraceHyperbolic,*)
(*but with an added termination condition \[CurlyPhi] at = \[Pi]/n*)
(*so that this need not be solved for later.*)


With[{zeta = \[FormalZeta]},
  zetaTraceCandidate[n_][a_] :=
    Module[{rhoClearance},
      rhoClearance = 10^-8;
      NDSolveValue[
        {
          zeta'[s] == zetaVel[n][a] @ zeta[s],
          zeta[0] == rhoSharp[n][a][0] + rhoClearance,
          WhenEvent[
            Or[
              Arg @ Conjugate @ zeta[s] > Pi / n,
              vi[n][a] @ zeta[s] < 0,
              Abs @ zeta[s] > 1
            ],
            "StopIntegration"
          ]
        }, zeta, {s, -Pi, 0},
        NoExtrapolation
      ]
    ];
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


poly[n_Integer] := RegularPolygon[{1, 0}, n];


polyVertices[n_Integer] := CirclePoints[{1, 0}, n];


(* ::Subsubsection:: *)
(*Complement of regular n-gon (z-space)*)


(* ::Text:: *)
(*This is quicker to plot than using RegionDifference with FullRegion.*)


polyComplement[n_Integer] :=
  Module[{mid},
    Table[
      (* Midpoint of k-th edge *)
      mid = Way[
        AngleVector[2 Pi (k - 1) / n],
        AngleVector[2 Pi k / n]
      ];
      (* Corresponding half space *)
      HalfSpace[-mid, mid]
    , {k, n}]
  ];


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


contStyle = LightGray;
streamStyle = LightGray;
nonStyle = Directive[Opacity[0.7], LightGray];
termStyle = Pink;
unphysStyle = Black;


ellStyle = Darker[Green];
hypStyle = Magenta;


upperStyle = Blue;
lowerStyle = Red;
convexStyle = Black;


pointStyle = PointSize[Large];
glowStyle = Directive[Thick, Yellow, Opacity[0.7]];


(* ::Section:: *)
(*Schwarz--Christoffel transformation*)


(* ::Subsection:: *)
(*Forward transformation: disk (\[Zeta]-space) to polygon (z-space)*)


(* ::Text:: *)
(*See https://mathematica.stackexchange.com/a/159086*)


Table[
  Module[
   {rMax,
    rNum, rValues,
    phiNum, phiValues
   },
    rMax = 1;
    Show[
      EmptyFrame[{-rMax, rMax}, {-rMax, rMax}],
      (* Equipotentials (T == const) *)
      rNum = 4;
      rValues = Subdivide[0, 1, rNum] // Most;
      ParametricPlot[
        Table[
          r Exp[I phi] // zMap[n] // ReIm
        , {r, rValues}] // Evaluate,
        {phi, 0, 2 Pi},
        PlotStyle -> contStyle
      ],
      (* Streamlines (parallel to \[Del]T) *)
      phiNum = 4;
      phiValues = Subdivide[0, 2 Pi, n phiNum] // Most;
      ParametricPlot[
        Table[
          r Exp[I phi] // zMap[n] // ReIm
        , {phi, phiValues}] // Evaluate,
        {r, 0, 1},
        PlotStyle -> streamStyle
      ],
      (* Unphysical domain *)
      Graphics @ {unphysStyle,
        polyComplement[n]
      }
    ]
  ] // Ex @ StringJoin["schwarz-christoffel-forward-", ToString[n], ".pdf"]
, {n, 3, 7}]


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
  , {n, 3, 5}]
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
, {n, 3, 5}]


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
, {n, 3, 5}]


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
, {n, 3, 5}]


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
, {n, 3, 5}]


(* ::Subsection:: *)
(*Small A plot (\[Zeta]-space)*)


(* ::Text:: *)
(*Plots to exaggerate \[CurlyPhi]-dependence.*)
(*We see that hyperbolic critical terminal points occur along \[CurlyPhi] = 2\[Pi] k/n,*)
(*whereas elliptic critical terminal points occur along \[CurlyPhi] = 2\[Pi] (k - 1/2)/n.*)


Table[
  Module[
   {a,
    eps, rhoN,
    rhoMax, rhoMaxUnphys,
    rhoMaxNon, rhoMinNon
   },
    (* Small value of A *)
    a = aNat[n][0] / 10^3;
    (* Plot *)
    eps = 0.1;
    rhoMax = 1;
    rhoMaxUnphys = 1 + eps;
    rhoMaxNon = rhoSharp[n][a][0];
    rhoMinNon = rhoSharp[n][a][Pi / n];
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
      ],
      (* Hyperbolic critical terminal points *)      
      (* (\[CurlyPhi] == 2 Pi k / n) *)
      Graphics @ {Directive[hypStyle, pointStyle],
        Point @ CirclePoints[{rhoMaxNon, 0}, n]
      },
      (* (corresponding contour) *)
      Graphics @ {hypStyle,
        Circle[{0, 0}, rhoMaxNon]
      },
      (* Elliptic critical terminal points *)
      (* (\[CurlyPhi] == 2 Pi (k - 1/2) / n) *)
      Graphics @ {Directive[ellStyle, pointStyle],
        Point @ CirclePoints[{rhoMinNon, Pi / n}, n]
      },
      (* (corresponding contour) *)
      Graphics @ {ellStyle,
        Circle[{0, 0}, rhoMinNon]
      }
    ]
  ] // Ex @ StringJoin["polygon_zeta-viable-small-a-", ToString[n],".pdf"]
, {n, 3, 5}]


(* ::Section:: *)
(*Traced boundary plots*)


(* ::Subsection:: *)
(*Hot regime*)


(* ::Subsubsection:: *)
(*Starting points (\[Zeta]-space)*)


Table[
  Module[
   {a, idList,
    eps, rhoMax, rhoMaxUnphys, rhoMaxNon
   },
    a = aHot[n];
    (* Group names *)
    idList = {"general", "hyperbolic"};
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
      (* Non-viable domain *)
      RegionPlot[vi[n][a][reZeta + I imZeta] < 0,
        {reZeta, -rhoMaxNon, rhoMaxNon},
        {imZeta, -rhoMaxNon, rhoMaxNon},
        BoundaryStyle -> termStyle,
        PlotStyle -> nonStyle
      ],
      (* Unphysical domain *)
      RegionPlot[RPolar[reZeta, imZeta] > 1,
        {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
        {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ],
      (* Starting points *)
      ListPlot[
        Table[startZetaHot[n][id] // ReIm, {id, idList}],
        LabelingFunction -> Function @ Placed[#2[[2]], Center],
        PlotLegends -> idList,
        PlotStyle -> Directive[Opacity[0.7], pointStyle]
      ]
    ]
  ] // Ex @ StringJoin["polygon_zeta-hot-traced-starting-", ToString[n],".pdf"]
, {n, 3, 5}]


(* ::Subsubsection:: *)
(*Starting points (z-space)*)


Table[
  Module[
   {idList,
    rMax,
    rNum, rValues,
    phiNum, phiValues
   },
    (* Group names *)
    idList = {"general", "hyperbolic"};
    (* Plot *)
    rMax = 1;
    Show[
      EmptyFrame[{-rMax, rMax}, {-rMax, rMax}],
      (* Equipotentials (T == const) *)
      rNum = 4;
      rValues = Subdivide[0, 1, rNum] // Most;
      ParametricPlot[
        Table[
          r Exp[I phi] // zMap[n] // ReIm
        , {r, rValues}] // Evaluate,
        {phi, 0, 2 Pi},
        PlotStyle -> contStyle
      ],
      (* Streamlines (parallel to \[Del]T) *)
      phiNum = 4;
      phiValues = Subdivide[0, 2 Pi, n phiNum] // Most;
      ParametricPlot[
        Table[
          r Exp[I phi] // zMap[n] // ReIm
        , {phi, phiValues}] // Evaluate,
        {r, 0, 1},
        PlotStyle -> streamStyle
      ],
      (* Unphysical domain *)
      Graphics @ {unphysStyle,
        polyComplement[n]
      },
      (* Starting points *)
      ListPlot[
        Table[startZetaHot[n][id] // zMap[n] // ReIm, {id, idList}],
        LabelingFunction -> Function @ Placed[#2[[2]], Center],
        PlotLegends -> idList,
        PlotStyle -> Directive[Opacity[0.7], pointStyle]
      ]
    ]
  ] // Ex @ StringJoin["polygon_z-hot-traced-starting-", ToString[n],".pdf"]
, {n, 3, 5}]


(* ::Subsubsection:: *)
(*General (\[Zeta]-space)*)


Table[
  Module[
   {a,
    eps, rhoMax, rhoMaxUnphys, rhoMaxNon
   },
    a = aHot[n];
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
        , {zeta, zetaTraHot[n][id]}]
      , {id, {"general"}}],
      Table[
        Table[
          Table[
            ParametricPlot[
              zeta[s] Exp[I 2 Pi k / n]
                // {#, Conjugate[#]} &
                // ReIm
                // Evaluate,
              {s, DomainStart[zeta], DomainEnd[zeta]},
              PlotStyle -> glowStyle
            ]
          , {k, 0, n - 1}]
        , {zeta, zetaTraHot[n][id]}]
      , {id, {"hyperbolic"}}],
      (* Unphysical domain *)
      RegionPlot[RPolar[reZeta, imZeta] > 1,
        {reZeta, -rhoMaxUnphys, rhoMaxUnphys},
        {imZeta, -rhoMaxUnphys, rhoMaxUnphys},
        BoundaryStyle -> None,
        PlotStyle -> unphysStyle
      ]
    ]
  ] // Ex @ StringJoin["polygon_zeta-hot-traced-full-", ToString[n],".pdf"]
, {n, 3, 5}]


(* ::Subsubsection:: *)
(*General (z-space)*)


Table[
  Module[
   {rMax,
    rNum, rValues,
    phiNum, phiValues
   },
    (* Plot *)
    rMax = 1;
    Show[
      EmptyFrame[{-rMax, rMax}, {-rMax, rMax}],
      (* Equipotentials (T == const) *)
      rNum = 4;
      rValues = Subdivide[0, 1, rNum] // Most;
      ParametricPlot[
        Table[
          r Exp[I phi] // zMap[n] // ReIm
        , {r, rValues}] // Evaluate,
        {phi, 0, 2 Pi},
        PlotStyle -> contStyle
      ],
      (* Streamlines (parallel to \[Del]T) *)
      phiNum = 4;
      phiValues = Subdivide[0, 2 Pi, n phiNum] // Most;
      ParametricPlot[
        Table[
          r Exp[I phi] // zMap[n] // ReIm
        , {phi, phiValues}] // Evaluate,
        {r, 0, 1},
        PlotStyle -> streamStyle
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
        , {zeta, zetaTraHot[n][id]}]
      , {id, {"general"}}],
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
              PlotStyle -> glowStyle
            ]
          , {k, 0, n - 1}]
        , {zeta, zetaTraHot[n][id]}]
      , {id, {"hyperbolic"}}],
      (* Unphysical domain *)
      Graphics @ {unphysStyle,
        polyComplement[n]
      }
    ]
  ] // Ex @ StringJoin["polygon_z-hot-traced-full-", ToString[n],".pdf"]
, {n, 3, 5}]


(* ::Subsubsection:: *)
(*Hyperbolic-only animations (\[Zeta]-space)*)


Table[
  Module[
   {ph,
    aN, aMin, aMax, aStep, aValues,
    eps, rhoMax, rhoMaxUnphys, rhoMaxNon,
    zeta
   },
    (* Azimuthal angle in \[Zeta]-space *)
    ph = 0;
    (* Values of A *)
    aN = aNat[n][ph];
    aMin = 1/10;
    aStep = 1/10;
    aMax = Floor[aN, aStep];
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
        ],
        (* Traced boundaries *)
        If[a < aN,
          (* If still in hot regime, plot *)
          zeta = zetaTraceHyperbolic[n][a];
          Table[
            ParametricPlot[
              zeta[s] Exp[I 2 Pi k / n]
                // {#, Conjugate[#]} &
                // ReIm
                // Evaluate,
              {s, DomainStart[zeta], DomainEnd[zeta]},
              PlotStyle -> hypStyle
            ]
          , {k, 0, n - 1}],
          (* Otherwise don't plot *)
          {}
        ]
      ]
    , {a, aValues}]
  ] // Ex[
    StringJoin["polygon_zeta-traced-hyperbolic-", ToString[n],".gif"],
    gifOpts
  ]
, {n, 3, 5}]


(* ::Subsubsection:: *)
(*Hyperbolic-only animations (z-space)*)


Table[
  Module[
   {ph,
    aN, aMin, aMax, aStep, aValues,
    rMax,
    rNum, rValues,
    phiNum, phiValues,
    zeta
   },
    (* Azimuthal angle in \[Zeta]-space *)
    ph = 0;
    (* Values of A *)
    aN = aNat[n][ph];
    aMin = 1/10;
    aStep = 1/10;
    aMax = Floor[aN, aStep];
    aValues = Range[aMin, aMax, aStep];
    aValues = Append[aValues, aN] // Sort[#, Less] &;
    (* Animation *)
    Table[
      rMax = 1;
      Show[
        EmptyFrame[{-rMax, rMax}, {-rMax, rMax},
          ImageSize -> 360,
          PlotLabel -> BoxedLabel[aIt == N[a]]
        ],
        (* Equipotentials (T == const) *)
        rNum = 4;
        rValues = Subdivide[0, 1, rNum] // Most;
        ParametricPlot[
          Table[
            r Exp[I phi] // zMap[n] // ReIm
          , {r, rValues}] // Evaluate,
          {phi, 0, 2 Pi},
          PlotStyle -> contStyle
        ],
        (* Streamlines (parallel to \[Del]T) *)
        phiNum = 4;
        phiValues = Subdivide[0, 2 Pi, n phiNum] // Most;
        ParametricPlot[
          Table[
            r Exp[I phi] // zMap[n] // ReIm
          , {phi, phiValues}] // Evaluate,
          {r, 0, 1},
          PlotStyle -> streamStyle
        ],
        (* Traced boundaries *)
        If[a < aN,
          (* If still in hot regime, plot *)
          zeta = zetaTraceHyperbolic[n][a];
          Table[
            ParametricPlot[
              zeta[s] Exp[I 2 Pi k / n]
                // zMap[n]
                // {#, Conjugate[#]} &
                // ReIm
                // Evaluate,
              {s, DomainStart[zeta], DomainEnd[zeta]},
              PlotStyle -> hypStyle
            ]
          , {k, 0, n - 1}],
          (* Otherwise don't plot *)
          {}
        ],
        (* Unphysical domain *)
        Graphics @ {unphysStyle,
          polyComplement[n]
        }
      ]
    , {a, aValues}]
  ] // Ex[
    StringJoin["polygon_z-traced-hyperbolic-", ToString[n],".gif"],
    gifOpts
  ]
, {n, 3, 5}]


(* ::Section:: *)
(*Traced boundary curvature*)
