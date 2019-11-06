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
  dest = "polygon-a-meeting.txt";
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
  zetaTraceCand[n_][a_] :=
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


(* ::Subsubsection:: *)
(*z' = dz/ds*)


(* ::Text:: *)
(*See (r4.39) (Page r4-11).*)


zVel[n_][a_][zeta_] :=
  Divide[
    I f[a][zeta] + Sqrt @ vi[n][a][zeta],
    gDer[n][zeta]
  ] // Evaluate;


(* ::Subsubsection:: *)
(*z'' = d^2(z)/ds^2*)


(* ::Text:: *)
(*See (r4.40) to (r4.43) (Pages r4-11 & r4-12).*)


zAcc[n_][a_][zeta_] := Plus[
  (* 1st term of (r4.40) *)
  Times[
    (* i F + sqrt(\[CapitalPhi]) *)
    I f[a][zeta] + Sqrt @ vi[n][a][zeta],
    (* d\[Zeta]/ds *)
    zetaVel[n][a][zeta],
    (* d/d\[Zeta] (1 / (dG/dz)) *)
    D[1 / gDer[n][zeta], zeta]
  ],
  (* 2nd term of (r4.40) *)
  Divide[
    Plus[
      (* i F' *)
      Times[
        I,
        -4 (Re @ g[zeta])^3 / a,
        Re @ (zetaVel[n][a][zeta] gDerZeta[zeta])
      ],
      (* (sqrt(\[CapitalPhi]))' \[Equal] \[CapitalPhi]' / (2 sqrt(\[CapitalPhi]))  *)
      Divide[
        Plus[
          (* 1st term of (r4.43) *)
          Times[
            2,
            Re @ gDer[n][zeta],
            Re @ (zetaVel[n][a][zeta] D[gDer[n][zeta], zeta])
          ],
          (* 2nd term of (r4.43) *)
          Times[
            2,
            Im @ gDer[n][zeta],
            Im @ (zetaVel[n][a][zeta] D[gDer[n][zeta], zeta])
          ],
          (* 3rd term of (r4.43) *)
          Times[
            -8 (Re @ g[zeta])^7 / a^2,
            Re @ (zetaVel[n][a][zeta] gDerZeta[zeta])
          ]
        ],
        2 Sqrt @ vi[n][a][zeta]
      ]
    ],
    (* dG/dz *)
    gDer[n][zeta]
  ]
] // Evaluate;


(* ::Subsubsection:: *)
(*Curvature x' y'' - y' x''*)


(* ::Text:: *)
(*See (r4.39a) (Page r4-11).*)


curTra[n_][a_][zeta_] :=
  Module[{xVel, yVel, xAcc, yAcc},
    {xVel, yVel} = zVel[n][a][zeta] // ReIm;
    {xAcc, yAcc} = zAcc[n][a][zeta] // ReIm;
    xVel yAcc - yVel xAcc
  ] // Evaluate;


(* ::Subsubsection:: *)
(*Curvature at \[CurlyPhi] = \[Pi]/n*)


(* ::Text:: *)
(*The curvature evaluated at the corners formed by the candidate boundaries.*)


curCandCorner[n_][a_] /; aMeet[n] <= a < aNat[n][0] :=
  Module[{zeta},
    zeta = zetaTraceCand[n][a];
    zeta @ DomainStart[zeta] // curTra[n][a]
  ];


(* ::Subsubsection:: *)
(*A_i (inflection dimensionless group)*)


(* ::Text:: *)
(*This is the smallest A > 0 for which \[CurlyPhi] = \[Pi]/n at \[Rho] = 1.*)
(*See Page r4-10.*)


(* Compute A_i using the bisection algorithm *)
aInfl = Module[
 {dest, nValues,
  aMin, aMax,
  a, num,
  aAss, numAss
 },
  (* (This is not slow, nevertheless compute once and store.) *)
  (* (Delete the file manually to compute from scratch.) *)
  dest = "polygon-a-inflection.txt";
  nValues = Range[3, 5];
  If[FileExistsQ[dest],
    (* If already stored, import *)
    {aAss, numAss} = Import[dest] // Uncompress,
    (* Otherwise compute and export for next time *)
    aAss = Association[];
    numAss = Association[];
    aMin = 15/10;
    aMax = 25/10;
    (* NOTE: these may not work for n > 5 *)
    Table[
      {a, num} = SeekRootBisection[
        curCandCorner[n][#] &,
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
  (* Return A_i values *)
  aAss
];
aInfl // N


(* ::Subsection:: *)
(*Representative values of A for a convex domain*)


aConvex[n_] := Ceiling[aInfl[n] + 0.25, 0.1];
Table[
  aConvex[n] < aNat[n][0]
, {n, 3, 5}]


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


(* ::Subsection:: *)
(*Repeated plots*)


(* ::Subsubsection:: *)
(*Equipotentials and streamlines (z-space)*)


equipStream[n_, opts : OptionsPattern[RegionPlot]] :=
  Module[{rMax, rhoNum, rhoValues, phNum, phValues},
    rMax = 1;
    Show[
      EmptyFrame[{-rMax, rMax}, {-rMax, rMax}, opts],
      (* Equipotentials (T == const) *)
      rhoNum = 4;
      rhoValues = Subdivide[0, 1, rhoNum] // Most;
      ParametricPlot[
        Table[
          rho Exp[I ph] // zMap[n] // ReIm
        , {rho, rhoValues}] // Evaluate,
        {ph, 0, 2 Pi},
        PlotStyle -> contStyle
      ],
      (* Streamlines (parallel to \[Del]T) *)
      phNum = 4;
      phValues = Subdivide[0, 2 Pi, n phNum] // Most;
      ParametricPlot[
        Table[
          rho Exp[I ph] // zMap[n] // ReIm
        , {ph, phValues}] // Evaluate,
        {rho, 0, 1},
        PlotStyle -> streamStyle
      ]
    ]
  ];


(* ::Subsubsection:: *)
(*Unphysical domain (z-space)*)


unphysDomain[n_] := Graphics @ {unphysStyle, polyComplement[n]};


(* ::Section:: *)
(*Schwarz--Christoffel transformation*)


(* ::Subsection:: *)
(*Forward transformation: disk (\[Zeta]-space) to polygon (z-space)*)


(* ::Text:: *)
(*See https://mathematica.stackexchange.com/a/159086*)


Table[
  Show[
    (* Equipotentials and streamlines *)
    equipStream[n],
    (* Unphysical domain *)
    unphysDomain[n]
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
  Module[{idList},
    (* Group names *)
    idList = {"general", "hyperbolic"};
    (* Plot *)
    Show[
      (* Equipotentials and streamlines *)
      equipStream[n],
      (* Unphysical domain *)
      unphysDomain[n],
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
  Module[{},
    Show[
      (* Equipotentials and streamlines *)
      equipStream[n],
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
      unphysDomain[n]
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
      Show[
        (* Equipotentials and streamlines *)
        equipStream[n,
          ImageSize -> 360,
          PlotLabel -> BoxedLabel[aIt == N[a]]
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
        unphysDomain[n]
      ]
    , {a, aValues}]
  ] // Ex[
    StringJoin["polygon_z-traced-hyperbolic-", ToString[n],".gif"],
    gifOpts
  ]
, {n, 3, 5}]


(* ::Subsection:: *)
(*Convex domain*)


(* ::Subsubsection:: *)
(*Without known solution (z-space)*)


Table[
  Module[{a, zeta, rMax},
    a = aConvex[n];
    zeta = zetaTraceCand[n][a];
    rMax = zeta @ DomainStart[zeta] // zMap[n] // 1.2 Abs[#] &;
    Show[
      (* Equipotentials and streamlines *)
      equipStream[n,
        PlotLabel -> BoxedLabel[aIt == N[a]],
        PlotRange -> {{-rMax, rMax}, {-rMax, rMax}}
      ],
      (* Traced boundaries *)
      Table[
        ParametricPlot[
          zeta[s] Exp[I 2 Pi k / n]
            // zMap[n]
            // {#, Conjugate[#]} &
            // ReIm
            // Evaluate,
          {s, DomainStart[zeta], DomainEnd[zeta]},
          PlotStyle -> convexStyle
        ]
      , {k, 0, n - 1}]
    ]
  ] // Ex @ StringJoin["polygon_z-convex-", ToString[n], ".pdf"]
, {n, 3, 5}]


(* ::Section:: *)
(*Traced boundary curvature*)


(* ::Subsection:: *)
(*Candidate boundaries*)


(* ::Subsubsection:: *)
(*Meeting (A = A_m) (z-space)*)


Table[
  Module[{a, zeta},
    a = aMeet[n];
    zeta = zetaTraceCand[n][a];
    Show[
      (* Equipotentials and streamlines *)
      equipStream[n,
        PlotLabel -> BoxedLabel[aIt == N[a]]
      ],
      (* Candidate traced boundaries *)
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
      (* Unphysical domain *)
      unphysDomain[n]
    ]
  ] // Ex @ StringJoin["polygon_z-candidate-meeting-", ToString[n], ".pdf"]
, {n, 3, 5}]


(* ::Subsubsection:: *)
(*Curvature at \[CurlyPhi] = \[Pi]/n*)


Module[
 {nValues, dest, aCurTables,
  aMin, aMax, num, aValues,
  inflDotStyle = Directive[Red, Opacity[0.7], pointStyle]
 },
  nValues = Range[3, 5];
  (* Compute tables of (A, cur.) values *)
  (* (This is slightly slow (~4 sec), so compute once and store.) *)
  (* (Delete the file manually to compute from scratch.) *)
  dest = "polygon-candidate-curvature-corner.txt";
  If[FileExistsQ[dest],
    (* If already stored, import *)
    aCurTables = Import[dest] // Uncompress,
    (* Otherwise compute and store for next time *)
    aCurTables = Table[
      aMin = aMeet[n];
      aMax = Min[aNat[n][0], 3];
      num = 100;
      aValues = Subdivide[aMin, aMax, num] // Most;
      Table[
        {a, curCandCorner[n][a]}
      , {a, aValues}]
    , {n, nValues}];
    aCurTables // Compress // Ex[dest]
  ];
  (* Plot *)
  Show[
    (* A vs corner curvature *)
    ListPlot[
      aCurTables,
      AxesLabel -> {aIt, "Corner cur."},
      Joined -> True,
      PlotLegends -> LineLegend[nValues, LegendLabel -> nIt],
      PlotRange -> {{All, 3}, {All, 5}},
      PlotOptions[Axes] // Evaluate
    ],
    (* A_i (inflection) *)
    Graphics @ {inflDotStyle,
      Point @ Table[{aInfl[n], 0}, {n, nValues}]
    },
    Graphics @ Text[
      Table[{n, aInfl[n]}, {n, nValues}]
        // TableForm[#, TableHeadings -> {{}, {nIt, Subscript[aIt, "i"]}}] &
        // Style[#, 17] &,
      {aMeet[3] + 0.15, 3}
    ]
  ] // Ex["polygon-candidate-curvature-corner.pdf"]
]


(* ::Subsubsection:: *)
(*Inflection (A = A_i) (z-space)*)


Table[
  Module[{a, zeta, rMax},
    a = aInfl[n];
    zeta = zetaTraceCand[n][a];
    rMax = zeta @ DomainStart[zeta] // zMap[n] // 1.2 Abs[#] &;
    Show[
      (* Equipotentials and streamlines *)
      equipStream[n,
        PlotLabel -> BoxedLabel[aIt == N[a]],
        PlotRange -> {{-rMax, rMax}, {-rMax, rMax}}
      ],
      (* Candidate traced boundaries *)
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
      , {k, 0, n - 1}]
    ]
  ] // Ex @ StringJoin["polygon_z-candidate-inflection-", ToString[n], ".pdf"]
, {n, 3, 5}]
