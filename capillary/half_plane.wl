(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< NDSolve`FEM`
<< Conway`
<< LaplaceYoung`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "half_plane"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Contact angles*)


(* gpd stands for "gamma per degree". *)
gpdValues = {1, 5, 10, 15, 30, 45, 60, 75};


(* ::Subsection:: *)
(*Finite element mesh*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["half_plane-mesh.txt",
  Module[
   {xMax, yMax, ySpacing, yValues,
    bPointList, nB, mod,
    bMesh, mesh,
    prWet
   },
    (* Length scales *)
    xMax = 10;
    yMax = 1/2;
    ySpacing = 0.005;
    (* Boundary points *)
    yValues = UniformRange[yMax, -yMax, -ySpacing];
    bPointList = Join[
      (* Fine length scale along x == 0 *)
      Table[{0, y}, {y, yValues}],
      (* Default length scale elsewhere *)
      {{xMax, -yMax}, {xMax, yMax}}
    ];
    nB = Length[bPointList];
    mod[n_] := Mod[#, n, 1] &;
    (* Boundary element mesh *)
    bMesh = ToBoundaryMesh[
      "Coordinates" -> bPointList,
      "BoundaryElements" -> {
        LineElement[
          Table[{n, n + 1}, {n, nB}] // mod[nB]
        ]
      }
    ];
    (* Build mesh *)
    mesh = ToElementMesh[bMesh,
      "ImproveBoundaryPosition" -> True
    ];
    (* Predicate function for wetting boundary (x == 0) *)
    prWet = Function[{x, y}, x == 0 // Evaluate];
    {mesh, prWet}
      // Compress
  ]
]


(* ::Subsection:: *)
(*Solve PDE (Version 12 required)*)


(* (This is slow (~7 sec), so compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["half_plane-solution.txt",
  Module[{mesh, prWet, gamma},
    (* Import mesh *)
    {mesh, prWet} = Import["half_plane-mesh.txt"] // Uncompress;
    (* Solve Laplace--Young equation *)
    Association @ Table[
      gamma = gpd * Degree;
      gpd -> SolveLaplaceYoung[gamma, mesh, prWet]
    , {gpd, gpdValues}]
      // Compress
  ]
]


(* ::Subsection:: *)
(*Solve PDE (fixed-point iteration)*)


(* (This is extremely slow (~10 min), so compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["half_plane-solution-fixed_point.txt",
  Module[{mesh, prWet, gamma},
    (* Import mesh *)
    {mesh, prWet} = Import["half_plane-mesh.txt"] // Uncompress;
    (* Solve Laplace--Young equation *)
    Association @ Table[
      gamma = gpd * Degree;
      gpd -> SolveLaplaceYoungFixedPoint[gamma, mesh, prWet]
    , {gpd, gpdValues}]
      // Compress
  ]
]


(* ::Subsection:: *)
(*Italicise symbols*)


gIt = Style["\[Gamma]"];


(* ::Subsection:: *)
(*Global styles for plots*)


pointStyle = PointSize[Large];


(* ::Section:: *)
(*Finite element mesh*)


Module[{mesh, prWet, nElem},
  (* Import mesh *)
  {mesh, prWet} = Import["half_plane-mesh.txt"] // Uncompress;
  nElem = Length @ mesh[[2, 1, 1]];
  (* Plot *)
  Show[
    mesh["Wireframe"],
    PlotLabel -> FString @ "{nElem} mesh elements",
    ImageSize -> 1440,
    PlotOptions[Axes] // Evaluate
  ]
] // Ex["half_plane-mesh.pdf"]


(* ::Section:: *)
(*Fine length scale*)


(* ::Text:: *)
(*Test what fine length scale to use along x = 0*)
(*by benchmarking against the exact half-plane solution for gamma = 1 degree.*)


(* (This is slow (~20 sec).) *)
Module[
  {
    gamma,
    hTheory,
    xMax, yMax,
    mod,
    prWet,
    ySpacingValues,
    yValues, boundaryPointList, numPoints,
    boundaryMesh, mesh,
    tSol, hNonlinear, relErrorNonlinear,
    dummyForTrailingCommas
  },
  (* Contact angle *)
  gamma = 1 Degree;
  (* Theoretical height rise *)
  hTheory = HHalfPlane[gamma];
  (* Length scales *)
  xMax = 10;
  yMax = 1/2;
  (* Modulo *)
  mod[n_] := Mod[#, n, 1] &;
  (* Predicate function for wetting boundary (x == 0) *)
  prWet = Function[{x, y}, x == 0];
  (* Fine length scales to test *)
  ySpacingValues = {0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1};
  (* For each fine length scale *)
  Table[
    (* Boundary points *)
    yValues = UniformRange[yMax, -yMax, -ySpacing];
    boundaryPointList = Join[
      (* Fine length scale along x == 0 *)
      Table[{0, y}, {y, yValues}],
      (* Default length scale elsewhere *)
      {{xMax, -yMax}, {xMax, yMax}},
      (* Dummy for trailing commas *)
      {}
    ];
    numPoints = Length[boundaryPointList];
    (* Boundary element mesh *)
    boundaryMesh =
      ToBoundaryMesh[
        "Coordinates" -> boundaryPointList,
        "BoundaryElements" -> {
          LineElement[
            Table[{n, n + 1}, {n, numPoints}] // mod[numPoints]
          ]
        },
        {}
      ];
    (* Build mesh *)
    mesh =
      ToElementMesh[
        boundaryMesh
        , "ImproveBoundaryPosition" -> True
    ];
    (* Compute wall height *)
    tSol = SolveLaplaceYoung[gamma, mesh, prWet];
    hNonlinear = tSol[0, 0];
    (* Compare with theoretical height *)
    relErrorNonlinear = hNonlinear / hTheory - 1;
    (* Return table row *)
    {
      ySpacing,
      Length @ mesh[[2, 1, 1]],
      hNonlinear,
      relErrorNonlinear,
      Nothing
    }
    , {ySpacing, ySpacingValues}
  ]
    //
      TableForm[#
        , TableHeadings -> {
            None,
            {
              "Fine length",
              "Mesh elements",
              "Computed height",
              "rel. error",
              Nothing
            }
          }
      ] &
    //
      Column @ {
        {"Contact angle", gamma},
        {"Theoretical height", hTheory // N},
        #
      } &
    //
      Ex["half_plane-fine-length-scale.pdf"]
]


(* ::Section:: *)
(*Numerical solution*)


(* ::Subsection:: *)
(*Using built-in nonlinear solver*)


Module[{tSolAss, tSol, mesh, gamma},
  tSolAss = Import["half_plane-solution.txt"] // Uncompress;
  Table[
    tSol = tSolAss[gpd];
    mesh = tSol["ElementMesh"];
    gamma = gpd * Degree;
    Block[{x = \[FormalX], y = \[FormalY]},
    (* (Using With results in SetDelayed::wrsym and an empty plot) *)
      Plot3D[tSol[x, y], Element[{x, y}, mesh],
        AxesLabel -> Italicise /@ {"x", "y", "T"},
        PlotLabel -> Column[
          {
            "Numerical solution",
            gIt == gamma
          },
          Alignment -> Center
        ],
        PlotRange -> {0, Sqrt[2]},
        PlotOptions[Axes] // Evaluate
      ]
    ] // Ex @ FString @ "half_plane-solution-gpd-{gpd}.png"
  , {gpd, gpdValues}]
]


(* ::Subsection:: *)
(*Using fixed-point iteration*)


Module[{tSolAss, tSol, n, mesh, gamma},
  tSolAss = Import["half_plane-solution-fixed_point.txt"] // Uncompress;
  Table[
    {tSol, n} = tSolAss[gpd];
    mesh = tSol["ElementMesh"];
    gamma = gpd * Degree;
    Block[{x = \[FormalX], y = \[FormalY]},
    (* (Using With results in SetDelayed::wrsym and an empty plot) *)
      Plot3D[tSol[x, y], Element[{x, y}, mesh],
        AxesLabel -> Italicise /@ {"x", "y", "T"},
        PlotLabel -> Column[
          {
            "Numerical solution: {n} iterations" // FString,
            gIt == gamma
          },
          Alignment -> Center
        ],
        PlotRange -> {0, Sqrt[2]},
        PlotOptions[Axes] // Evaluate
      ]
    ] // Ex @ FString @ "half_plane-solution-fixed_point-gpd-{gpd}.png"
  , {gpd, gpdValues}]
]


(* ::Subsection:: *)
(*Discrepancy between nonlinear solver and fixed-point iteration*)


Module[
 {tSolAssNonlinear, tSolAssFixedPoint,
  tSolNonlinear, tSolFixedPoint,
  gamma, xMax
 },
  tSolAssNonlinear = Import["half_plane-solution.txt"] // Uncompress;
  tSolAssFixedPoint = Import["half_plane-solution-fixed_point.txt"] // Uncompress;
  Table[
    tSolNonlinear = tSolAssNonlinear[gpd];
    tSolFixedPoint = tSolAssFixedPoint[gpd] // First;
    gamma = gpd * Degree;
    xMax = 3;
    Plot[tSolNonlinear[x, 0] - tSolFixedPoint[x, 0], {x, 0, xMax},
      AxesLabel -> {Italicise @ "x", ""},
      PlotLabel -> Column[
        {
          "[Nonlinear]" - "[Fixed-point]",
          gIt == gamma
        },
        Alignment -> Center
      ],
      PlotRange -> Full,
      PlotStyle -> Red,
      PlotOptions[Axes] // Evaluate
    ] // Ex @ FString @ "half_plane-solution-discrepancy-gpd-{gpd}.pdf"
  , {gpd, gpdValues}]
]


(* ::Subsection:: *)
(*Compare with theoretical height rise values*)


Module[
 {tSolAssNonlinear, tSolAssFixedPoint,
  tSolNonlinear, tSolFixedPoint, n,
  gamma, hTheory,
  hNonlinear, relErrorNonlinear,
  hFixedPoint, relErrorFixedPoint
 },
  tSolAssNonlinear = Import["half_plane-solution.txt"] // Uncompress;
  tSolAssFixedPoint = Import["half_plane-solution-fixed_point.txt"] // Uncompress;
  Table[
    tSolNonlinear = tSolAssNonlinear[gpd];
    {tSolFixedPoint, n} = tSolAssFixedPoint[gpd];
    (* Theoretical *)
    gamma = gpd * Degree;
    hTheory = HHalfPlane[gamma];
    (* Nonlinear solver *)
    hNonlinear = tSolNonlinear[0, 0];
    relErrorNonlinear = hNonlinear / hTheory - 1;
    (* Fixed-point iteration *)
    hFixedPoint = tSolFixedPoint[0, 0];
    relErrorFixedPoint = hFixedPoint / hTheory - 1;
    {
      gamma, hTheory // N,
      hNonlinear, relErrorNonlinear,
      hFixedPoint, relErrorFixedPoint,
      Abs[relErrorNonlinear] < Abs[relErrorFixedPoint]
    }
  , {gpd, gpdValues}]
] // TableForm[#,
  TableHeadings -> {
    None,
    {
      "\[Gamma]", "Theory",
      "Nonlinear", "rel. error",
      "Fixed-point", "rel. error",
      "Nonlinear better"
    }
  }
] & // Ex["half_plane-height-comparison.pdf"]


(* ::Subsection:: *)
(*Numerical error for nonlinear solver (\[Gamma] = 1\[Degree])*)


Module[
 {gpd, gamma,
  tSolAss, tSol,
  xFun, xMax, tMax, tMin, tStep, tInterp
 },
  (* Contact angle *)
  gpd = 1;
  gamma = gpd * Degree;
  (* Import numerical solution *)
  tSolAss = Import["half_plane-solution.txt"] // Uncompress;
  tSol = tSolAss[gpd];
  (* Implicit exact half-plane solution x == x(T) *)
  xFun[t_] := XHalfPlane[gamma][t] // Evaluate;
  (* Interpolation to obtain T == T(x) *)
  tMax = HHalfPlane[gamma];
  xMax = 3;
  tMin = Quiet[
    1/2 * SeekRoot[xFun[#] - xMax &, {0, tMax}]
  , {Power::infy}];
  tStep = 0.001; (* Recall ySpacing = 0.005 above *)
  tInterp = Interpolation @ (
    Table[{xFun[t], t}, {t, tMax, tMin, -tStep}]
  );
  (* Plot *)
  LogPlot[tSol[x, 0] / tInterp[x] - 1 // Abs,
    {x, 0, xMax},
    AxesLabel -> {Italicise @ "x", "Rel. error"},
    AxesOrigin -> {-0.125, Automatic},
    Exclusions -> None,
    PlotRange -> Full,
    PlotOptions[Axes] // Evaluate
  ] // Ex @ FString["half_plane-solution-rel_error-gpd-{gpd}.pdf"]
]
