(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< NDSolve`FEM`
<< Conway`
<< Curvilinear`
<< LaplaceYoung`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "wedge_obtuse"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Finite element mesh*)


(* Wedge half-angles (apd means "alpha per degree") *)
apdMeshValues = {120, 135, 150};
rMaxMesh = 10;


(* (These are not slow, nevertheless compute once and store.) *)
(* (Delete the files manually to compute from scratch.) *)
Table[
  ExportIfNotExists[FString @ "mesh/wedge_obtuse-mesh-apd-{apd}.txt",
    Module[
      {
        rMax, alpha,
        fineLengthScale, coarseLengthScale,
        boundaryPoints, numBoundaryPoints, boundaryElements,
        boundaryMesh, mesh,
        predicateWet,
        dummyForTrailingCommas
      },
      (* Geometry *)
      rMax = rMaxMesh;
      alpha = apd * Degree;
      (* Refinement length scales *)
      fineLengthScale = 0.01;
      coarseLengthScale = 0.2;
      (* Boundary points and elements *)
      boundaryPoints =
        Join[
          (* Lower wedge wall *)
          Table[
            XYPolar[r, -alpha]
            , {r, UniformRange[0, rMax, fineLengthScale]}
          ]
            // Most,
          (* Far arc at infinity *)
          Table[
            XYPolar[rMax, phi]
            , {phi, UniformRange[-alpha, alpha, coarseLengthScale / rMax]}
          ]
            // Most,
          (* Upper wedge wall *)
          Table[
            XYPolar[r, +alpha]
            , {r, UniformRange[rMax, 0, -fineLengthScale]}
          ]
            // Most,
          {}
        ];
      numBoundaryPoints = Length[boundaryPoints];
      boundaryElements =
        LineElement /@ {
          Table[{n, n + 1}, {n, numBoundaryPoints}]
            // Mod[#, numBoundaryPoints, 1] &
        };
      (* Build mesh *)
      boundaryMesh =
        ToBoundaryMesh[
          "Coordinates" -> boundaryPoints,
          "BoundaryElements" -> boundaryElements,
          {}
        ];
      mesh = ToElementMesh[boundaryMesh
        , "ImproveBoundaryPosition" -> True
        , MeshRefinementFunction -> MeshRefinementUniform[coarseLengthScale]
      ];
      (* Predicate function for wetting boundaries *)
      predicateWet =
        Function[{x, y},
          Abs[y] <= rMax Sin[Pi - alpha] && x <= 0
          // Evaluate
        ];
      (* Return *)
      {mesh, predicateWet} // Compress
    ]
  ]
  , {apd, apdMeshValues}
]


(* ::Subsection:: *)
(*Solve capillary BVP*)


(* Solution contact angles (gpd means "gamma per degree") *)
gpdSolutionValues = Join[{1}, Range[5, 85, 5]];


(* (These are slow (~10 min), so compute once and store.) *)
(* (Delete the files manually to compute from scratch.) *)
(*
  I do not know how much memory is required for this,
  so I would suggest running this on a fresh kernel.
*)
(* For each value of alpha: *)
Table[
  (* If all values of contact angle have been account for: *)
  If[
    Table[
      FString @ "solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"
      , {gpd, gpdSolutionValues}
    ] // AllTrue[FileExistsQ]
    ,
    (* Then skip: *)
    Null
    ,
    (* Otherwise: *)
    Module[{mesh, predicateWet, gamma, evaluationData},
      (* Import mesh *)
      {mesh, predicateWet} =
        FString @ "mesh/wedge_obtuse-mesh-apd-{apd}.txt"
          // Import // Uncompress;
      (* For each value of gamma: *)
      Table[
        gamma = gpd * Degree;
        ExportIfNotExists[
          FString @ "solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"
          ,
          evaluationData = EvaluationData @ SolveLaplaceYoung[gamma, mesh, predicateWet];
          evaluationData /@ {"Result", "Success", "FailureType", "MessagesText", "AbsoluteTiming"}
            // Compress
        ]
        , {gpd, gpdSolutionValues}
      ]
    ]
  ]
  , {apd, apdMeshValues}
]


(* ::Section:: *)
(*Finite element mesh check*)


(* (This is slow (~5 sec).) *)
Table[
  Module[
    {
      mesh, predicateWet, numElements,
      boundaryCoordinates, wetBoundaryCoordinates,
      plot,
      dummyForTrailingCommas
    },
    (* Import mesh *)
    {mesh, predicateWet} =
      FString @ "mesh/wedge_obtuse-mesh-apd-{apd}.txt"
        // Import // Uncompress;
    numElements = Length @ mesh[[2, 1, 1]];
    (* Determine predicate-satisfying boundary coordinates *)
    boundaryCoordinates =
      Part[
        mesh["Coordinates"],
        mesh["BoundaryElements"][[1, 1]] // Flatten // Union
      ];
    wetBoundaryCoordinates =
      Select[
        boundaryCoordinates,
        predicateWet @@ # &
      ];
    (* Make plot *)
    plot =
      Show[
        mesh["Wireframe"],
        ListPlot[boundaryCoordinates, PlotStyle -> Directive[Red, PointSize[Medium]]],
        ListPlot[wetBoundaryCoordinates, PlotStyle -> Directive[Blue, PointSize[Medium]]],
        {}
        , PlotLabel -> FString @ "{numElements} mesh elements"
      ];
    (* Export *)
    {
      (* Full plot *)
      plot
        // Ex @ FString @ "mesh/wedge_obtuse-mesh-apd-{apd}.png"
      ,
      (* Zoom near origin to verify wetting predicate *)
      Show[plot
        , ImageSize -> 480
        , PlotLabel -> None
        , PlotRange -> {{-0.15, 0.15}, {-0.15, 0.15}}
      ]
        // Ex @ FString @ "mesh/wedge_obtuse-mesh-apd-{apd}-zoom.png"
    }
  ]
  , {apd, apdMeshValues}
]


(* ::Section:: *)
(*Numerical solution*)


(* ::Subsection:: *)
(*Table*)


(* (This is slow (~30 sec) from all the calls to Import.) *)
Join @@ Table[
  Table[
    Module[{tNumerical, messages, timing, heightCorner},
      (* Import numerical solution etc. *)
      {tNumerical, messages, timing} =
        FString @ "solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"
          // Import // Uncompress // Part[#, {1, 4, 5}] &;
      (* Filter out uninformative messages *)
      messages = DeleteCases[messages, string_ /; StringContainsQ[string, "StringForm"]];
      (* Make pretty *)
      messages = StringRiffle[messages, "\n"];
      messages = messages /. {"" -> "\[HappySmiley]"};
      (* Compute corner height *)
      heightCorner = tNumerical[0, 0];
      (* Return row of table *)
      {apd * Degree, gpd * Degree, heightCorner, timing, messages}
    ]
    , {gpd, gpdSolutionValues}
  ]
  , {apd, apdMeshValues}
] // TableForm[#
  , TableAlignments -> {Left, Center}
  , TableDepth -> 2
  , TableHeadings -> {
      None,
      {"alpha", "gamma", "corner height", "absolute time", "messages"}
    }
] & // Ex["wedge_obtuse-solution-table.pdf"]


(* ::Section:: *)
(*Slope dodginess test*)


(*
  Numerical slope near a re-entrant corner is dodgy.
  Here, we examine the effect of additional refinement for:
  * alpha == {45, 135} Degree
  * gamma == 60 Degree
  * global coarse length scale 0.2
  * wall fine length scale 0.01
  * near-corner (r < 0.01) fine length scales {0.01, 0.005, 0.002, 0.001}
*)


(* ::Subsection:: *)
(*Compute data and store*)


(* (This is slow (~1.5 min).) *)
(*
  I do not know how much memory is required for this,
  so I would suggest running this on a fresh kernel.
*)
Module[
  {
    gamma,
    globalCoarseLengthScale, wallFineLengthScale, nearCornerRadius,
    apdValues, nearCornerFineLengthScaleValues,
      rMax, alpha,
      boundaryPoints, numBoundaryPoints, boundaryElements,
      equilateralTriangleArea, refinementFunction,
      boundaryMesh, mesh, numMeshElements,
      predicateWet,
      tNumerical,
      xValues, heightData, slopeData,
    dummyForTrailingCommas
  },
  (* Constants *)
  gamma = 60 Degree;
  globalCoarseLengthScale = 0.2;
  wallFineLengthScale = 0.01;
  nearCornerRadius = 0.01;
  (* Values to test *)
  apdValues = {45, 135};
  nearCornerFineLengthScaleValues = {0.01, 0.001, 0.0005, 0.0002, 0.0001};
  (* For each pair of test values: *)
  Table[
    (* Build mesh *)
    rMax = rMaxMesh;
    alpha = apd * Degree;
    (* Boundary points and elements *)
    boundaryPoints =
      Join[
        (* Lower wedge wall *)
        Table[
          XYPolar[r, -alpha]
          , {r, UniformRange[0, rMax, wallFineLengthScale]}
        ]
          // Most,
        (* Far arc at infinity *)
        Table[
          XYPolar[rMax, phi]
          , {phi, UniformRange[-alpha, alpha, globalCoarseLengthScale / rMax]}
        ]
          // Most,
        (* Upper wedge wall *)
        Table[
          XYPolar[r, +alpha]
          , {r, UniformRange[rMax, 0, -wallFineLengthScale]}
        ]
          // Most,
        {}
      ];
    numBoundaryPoints = Length[boundaryPoints];
    boundaryElements =
      LineElement /@ {
        Table[{n, n + 1}, {n, numBoundaryPoints}]
          // Mod[#, numBoundaryPoints, 1] &
      };
    (* Special refinement *)
    equilateralTriangleArea = Function[{len}, Sqrt[3]/4 * len^2];
    refinementFunction =
      Function[{vertices, area},
        Or[
          area > equilateralTriangleArea[globalCoarseLengthScale],
          Norm @ Mean[vertices] < nearCornerRadius
            && area > equilateralTriangleArea[nearCornerFineLengthScale]
        ]
      ];
    (* Build mesh *)
    boundaryMesh =
      ToBoundaryMesh[
        "Coordinates" -> boundaryPoints,
        "BoundaryElements" -> boundaryElements,
        {}
      ];
    mesh = ToElementMesh[boundaryMesh
      , "ImproveBoundaryPosition" -> True
      , MeshRefinementFunction -> refinementFunction
    ];
    (* Export mesh wireframe for visual check *)
    numMeshElements = Length @ mesh[[2, 1, 1]];
    Show[
      mesh["Wireframe"]
      , Frame -> True
      , ImageSize -> 240
      , PlotLabel -> Column @ {
          FString @ "apd: {apd}",
          FString @ "nearCornerFineLengthScale: {nearCornerFineLengthScale}",
          FString @ "numMeshElements: {numMeshElements}",
          Nothing
        }
      , PlotRange -> 0.1 {{-1, 1}, {-1, 1}}
      , PlotRangeClipping -> True
    ] // Ex @ FString[
      "slope_test/wedge_obtuse-slope_test-apd-{apd}-fls-{nearCornerFineLengthScale}.png"
    ];
    (* Predicate function for wetting boundaries *)
    predicateWet =
      If[alpha < Pi/2,
        Function[{x, y}, x <= rMax Cos[alpha] // Evaluate],
        Function[{x, y}, Abs[y] <= rMax Sin[Pi - alpha] && x <= 0 // Evaluate]
      ];
    (* Compute numerical solution to capillary BVP *)
    tNumerical = SolveLaplaceYoung[gamma, mesh, predicateWet];
    (* Compute data *)
    xValues =
      Join[
        Range[0, 0.1, 0.001],
        Range[0.1, 1, 0.01] // Rest,
        Range[1, 5, 0.1] // Rest,
        {}
      ];
    heightData = Table[{x, tNumerical[x, 0]}, {x, xValues}];
    slopeData = Table[{x, Derivative[1, 0][tNumerical][x, 0]}, {x, xValues}];
    (* Return row *)
    {apd, nearCornerFineLengthScale, heightData, slopeData}
    , {apd, apdValues}
    , {nearCornerFineLengthScale, nearCornerFineLengthScaleValues}
  ] // Compress
] // Ex["slope_test/wedge_obtuse-slope_test-all-data.txt"]
