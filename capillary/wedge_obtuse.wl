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


(* ::Subsection:: *)
(*Elliptic critical terminal point*)


x0[tNumerical_, gammaTracing_] :=
  Module[{p},
    (* \[PartialD]T/\[PartialD]x *)
    p = Derivative[1, 0][tNumerical];
    (* Find x == x_0 such that \[PartialD]T/\[PartialD]x == cot(gamma-tracing) *)
    SeekRoot[
      p[#, 0] + Cot[gammaTracing] &,
      {0, rMaxMesh}, 20
    ]
  ];


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


(* ::Subsection:: *)
(*Table for theoretical vs computed corner slopes*)


(* (This is slow (~10 sec) from all the calls to Import.) *)
Module[
  {
    apd, alpha,
    gamma, k,
    infiniteSlopeQ, theoreticalSlope,
    tNumerical, tXDerivative, tYDerivative,
    computedSlope, slopeRelError,
    dummyForTrailingCommas
  },
  (* Wedge half-angle *)
  apd = 135;
  alpha = apd * Degree;
  (* For each contact angle: *)
  Table[
    (* Theoretical corner slope *)
    gamma = gpd Degree;
    k = Sin[alpha] / Cos[gamma];
    infiniteSlopeQ = alpha >= Pi/2 + gamma;
    theoreticalSlope = If[infiniteSlopeQ, Infinity, 1 / Sqrt[k^2 - 1] // N];
    (* Import numerical solution *)
    tNumerical =
      Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
        // Uncompress // First;
    tXDerivative = Derivative[1, 0][tNumerical];
    tYDerivative = Derivative[0, 1][tNumerical];
    computedSlope = Sqrt[tXDerivative[0, 0]^2 + tYDerivative[0, 0]^2];
    slopeRelError = If[infiniteSlopeQ, "N/A", computedSlope / theoreticalSlope - 1];
    (* Build row of table *)
    {gamma, theoreticalSlope, computedSlope, slopeRelError}
    , {gpd, gpdSolutionValues}
  ]
    // TableForm[#
      , TableHeadings -> {
          None,
          {"gamma", "theoretical slope", "computed slope", "rel. error"}
        }
    ] &
] // Ex["wedge_obtuse-slope-comparison.pdf"]


(* ::Subsection:: *)
(*Corner neighbourhood test*)


(* A visual check of local linearity. *)
Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical, tXDerivative,
    coordMax, coordSliceValues,
    imageSize,
    dummyForTrailingCommas
  },
  (* Parameters *)
  apd = 135;
  gpd = 60;
  alpha = apd * Degree;
  gamma = gpd * Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  tXDerivative = Derivative[1, 0][tNumerical];
  (* Plot *)
  coordMax = 0.01;
  coordSliceValues = Subdivide[-coordMax, coordMax, 4];
  imageSize = 360;
  {
    (* Height *)
    Plot[
      Table[tNumerical[x, y], {y, coordSliceValues}] // Evaluate
      , {x, -coordMax, coordMax}
      , AxesLabel -> {"x", "height"}
      , AxesOrigin -> {-coordMax, Automatic}
      , ImageSize -> imageSize
      , PlotLegends -> LineLegend[coordSliceValues
          , LegendLabel -> Italicise["y"]
        ]
      , PlotRange -> Full
      , PlotStyle -> {Automatic, Automatic, Automatic, Dotted, Dotted}
      , PlotOptions[Axes] // Evaluate
    ]
      // Ex["solution/wedge_obtuse-solution-neighbourhood-height.png"]
    ,
    (* x-derivative *)
    Plot[
      Table[tXDerivative[x, y], {y, coordSliceValues}] // Evaluate
      , {x, -coordMax, coordMax}
      , AxesLabel -> {"x", "slope"}
      , AxesOrigin -> {-coordMax, Automatic}
      , ImageSize -> imageSize
      , PlotLegends -> LineLegend[coordSliceValues
          , LegendLabel -> Italicise["y"]
        ]
      , PlotRange -> Full
      , PlotStyle -> {Automatic, Automatic, Automatic, Dotted, Dotted}
      , PlotOptions[Axes] // Evaluate
    ]
      // Ex["solution/wedge_obtuse-solution-neighbourhood-slope.png"]
  }
]


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
      tNumerical, tXDerivative, tYDerivative, tSlope,
      rValues,
        symmetryHeightData, symmetrySlopeData,
        wallHeightData, wallSlopeData,
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
    tXDerivative = Derivative[1, 0][tNumerical];
    tYDerivative = Derivative[0, 1][tNumerical];
    tSlope = Function[{x, y}, Sqrt[tXDerivative[x, y]^2 + tYDerivative[x, y]^2]];
    rValues =
      Join[
        Range[0, 0.1, 0.001],
        Range[0.1, 1, 0.01] // Rest,
        Range[1, 5, 0.1] // Rest,
        {}
      ];
    (* Compute data (along line of symmetry) *)
    symmetryHeightData = Table[{r, tNumerical @@ XYPolar[r, 0]}, {r, rValues}];
    symmetrySlopeData = Table[{r, tSlope @@ XYPolar[r, 0]}, {r, rValues}];
    (* Compute data (along wall) *)
    wallHeightData = Table[{r, tNumerical @@ XYPolar[r, alpha]}, {r, rValues}];
    wallSlopeData = Table[{r, tSlope @@ XYPolar[r, alpha]}, {r, rValues}];
    (* Return row *)
    {
      apd, nearCornerFineLengthScale,
      symmetryHeightData, symmetrySlopeData,
      wallHeightData, wallSlopeData,
      Nothing
    }
    , {apd, apdValues}
    , {nearCornerFineLengthScale, nearCornerFineLengthScaleValues}
  ] // Compress
] // Ex["slope_test/wedge_obtuse-slope_test-all-data.txt"]


(* ::Subsection:: *)
(*Visualise*)


Module[
  {
    allData,
    apdValues, nearCornerFineLengthScaleValues,
    symmetryHeightData, symmetrySlopeData,
    cornerHeight, cornerSlope,
    wallHeightData, wallSlopeData,
    dummyForTrailingCommas
  },
  (* Import all data *)
  allData = Import["slope_test/wedge_obtuse-slope_test-all-data.txt"] // Uncompress;
  allData = Join @@ allData;
  (* Extract data *)
  apdValues = allData[[All, 1]] // Union;
  nearCornerFineLengthScaleValues = allData[[All, 2]] // Union;
  allData // Cases[
    {apd_, ncfls_, shData_, ssData_, whData_, wsData_} :> (
      (* Along line of symmetry *)
      symmetryHeightData[apd, ncfls] = shData;
        cornerHeight[apd, ncfls] = shData[[1, 2]];
      symmetrySlopeData[apd, ncfls] = ssData;
        cornerSlope[apd, ncfls] = ssData[[1, 2]];
      (* Along wall *)
      wallHeightData[apd, ncfls] = whData;
      wallSlopeData[apd, ncfls] = wsData;
    )
  ];
  (* For each wedge half-angle: *)
  Table[
    List[
      (* Corner height plot *)
      ListLogLinearPlot[
        Table[{ncfls, cornerHeight[apd, ncfls]}, {ncfls, nearCornerFineLengthScaleValues}]
        , AxesLabel -> {"ncfls", "corner height"}
        , ImageSize -> 360
        , Joined -> True
        , PlotLabel -> {"apd", apd}
        , PlotMarkers -> {Automatic, Medium}
        , PlotRange -> All
        , PlotOptions[Axes] // Evaluate
      ]
        // Ex @ FString["slope_test/wedge_obtuse-slope_test-corner-height-apd-{apd}.png"]
      ,
      (* Symmetry height plot *)
      ListPlot[
        Table[symmetryHeightData[apd, ncfls], {ncfls, nearCornerFineLengthScaleValues}]
        , AxesLabel -> {"r", "symm height"}
        , ImageSize -> 360
        , Joined -> True
        , PlotLabel -> {"apd", apd}
        , PlotLegends -> LineLegend[
            nearCornerFineLengthScaleValues
            , LegendLabel -> "ncfls"
            , LegendLayout -> "ReversedColumn"
          ]
        , PlotMarkers -> {Automatic, Medium}
        , PlotRange -> {{0, 0.01}, All}
        , PlotRangeClipping -> False
        , PlotOptions[Axes] // Evaluate
      ]
        // Ex @ FString["slope_test/wedge_obtuse-slope_test-symm-height-apd-{apd}.png"]
      ,
      (* Corner slope plot *)
      ListLogLinearPlot[
        Table[{ncfls, -cornerSlope[apd, ncfls]}, {ncfls, nearCornerFineLengthScaleValues}]
        , AxesLabel -> {"ncfls", "corner slope"}
        , ImageSize -> 360
        , Joined -> True
        , PlotLabel -> {"apd", apd}
        , PlotMarkers -> {Automatic, Medium}
        , PlotRange -> All
        , PlotOptions[Axes] // Evaluate
      ]
        // Ex @ FString["slope_test/wedge_obtuse-slope_test-corner-slope-apd-{apd}.png"]
      ,
      (* Symmetry slope plot *)
      ListPlot[
        Table[symmetrySlopeData[apd, ncfls], {ncfls, nearCornerFineLengthScaleValues}]
        , AxesLabel -> {"r", "symm slope"}
        , ImageSize -> 360
        , Joined -> True
        , PlotLabel -> {"apd", apd}
        , PlotLegends -> LineLegend[
            nearCornerFineLengthScaleValues
            , LegendLabel -> "ncfls"
            , LegendLayout -> "ReversedColumn"
          ]
        , PlotMarkers -> {Automatic, Medium}
        , PlotRange -> {{0, 0.01}, All}
        , PlotRangeClipping -> False
        , PlotRangePadding -> {None, Scaled[0.1]}
        , PlotOptions[Axes] // Evaluate
      ]
        // Ex @ FString["slope_test/wedge_obtuse-slope_test-symm-slope-apd-{apd}.png"]
      ,
      (* Wall height *)
      ListPlot[
        Table[wallHeightData[apd, ncfls], {ncfls, nearCornerFineLengthScaleValues}]
        , AxesLabel -> {"r", "wall height"}
        , ImageSize -> 360
        , Joined -> True
        , PlotLabel -> {"apd", apd}
        , PlotLegends -> LineLegend[
            nearCornerFineLengthScaleValues
            , LegendLabel -> "ncfls"
            , LegendLayout -> "ReversedColumn"
          ]
        , PlotMarkers -> {Automatic, Medium}
        , PlotRange -> {{0, 0.02}, All}
        , PlotRangeClipping -> False
        , PlotOptions[Axes] // Evaluate
      ]
        // Ex @ FString["slope_test/wedge_obtuse-slope_test-wall-height-apd-{apd}.png"]
      ,
      (* Wall slope *)
      ListPlot[
        Table[wallSlopeData[apd, ncfls], {ncfls, nearCornerFineLengthScaleValues}]
        , AxesLabel -> {"r", "wall slope"}
        , ImageSize -> 360
        , Joined -> True
        , PlotLabel -> {"apd", apd}
        , PlotLegends -> LineLegend[
            nearCornerFineLengthScaleValues
            , LegendLabel -> "ncfls"
            , LegendLayout -> "ReversedColumn"
          ]
        , PlotMarkers -> {Automatic, Medium}
        , PlotRange -> {{0, 0.02}, All}
        , PlotRangeClipping -> False
        , PlotOptions[Axes] // Evaluate
      ]
        // Ex @ FString["slope_test/wedge_obtuse-slope_test-wall-slope-{apd}.png"]
    ]
    , {apd, apdValues}
  ]
]


(* ::Subsection:: *)
(*Corner neighbourhood test*)


(* (This is slow (~15 sec).) *)
(*
  I do not know how much memory is required for this,
  so I would suggest running this on a fresh kernel.
*)
Module[
  {
    apd, gpd,
    alpha, gamma,
      globalCoarseLengthScale, wallFineLengthScale, nearCornerRadius,
      nearCornerFineLengthScale,
      rMax,
      boundaryPoints, numBoundaryPoints, boundaryElements,
      equilateralTriangleArea, refinementFunction,
      boundaryMesh, mesh, numMeshElements,
      predicateWet,
    tNumerical, tXDerivative, tYDerivative,
    coordMax, coordSliceValues,
    imageSize,
    dummyForTrailingCommas
  },
  (* Parameters *)
  apd = 135;
  gpd = 60;
  alpha = apd * Degree;
  gamma = gpd * Degree;
  (* Build mesh and compute numerical solution *)
  (* (copied from above) *)
    (* Parameters *)
    globalCoarseLengthScale = 0.2;
    wallFineLengthScale = 0.01;
    nearCornerRadius = 0.01;
    nearCornerFineLengthScale = 0.0001;
    (* Geometry *)
    rMax = rMaxMesh;
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
    (* Predicate function for wetting boundaries *)
    predicateWet =
      If[alpha < Pi/2,
        Function[{x, y}, x <= rMax Cos[alpha] // Evaluate],
        Function[{x, y}, Abs[y] <= rMax Sin[Pi - alpha] && x <= 0 // Evaluate]
      ];
    (* Compute numerical solution to capillary BVP *)
    tNumerical = SolveLaplaceYoung[gamma, mesh, predicateWet];
    tXDerivative = Derivative[1, 0][tNumerical];
  (* Plot *)
  coordMax = 0.01;
  coordSliceValues = Subdivide[-coordMax, coordMax, 4];
  imageSize = 360;
  {
    (* Height *)
    Plot[
      Table[tNumerical[x, y], {y, coordSliceValues}] // Evaluate
      , {x, -coordMax, coordMax}
      , AxesLabel -> {"x", "height"}
      , AxesOrigin -> {-coordMax, Automatic}
      , ImageSize -> imageSize
      , PlotLegends -> LineLegend[coordSliceValues
          , LegendLabel -> Italicise["y"]
        ]
      , PlotRange -> Full
      , PlotStyle -> {Automatic, Automatic, Automatic, Dotted, Dotted}
      , PlotOptions[Axes] // Evaluate
    ]
      // Ex["slope_test/wedge_obtuse-slope_test-neighbourhood-height.png"]
    ,
    (* x-derivative *)
    Plot[
      Table[tXDerivative[x, y], {y, coordSliceValues}] // Evaluate
      , {x, -coordMax, coordMax}
      , AxesLabel -> {"x", "slope"}
      , AxesOrigin -> {-coordMax, Automatic}
      , ImageSize -> imageSize
      , PlotLegends -> LineLegend[coordSliceValues
          , LegendLabel -> Italicise["y"]
        ]
      , PlotRange -> Full
      , PlotStyle -> {Automatic, Automatic, Automatic, Dotted, Dotted}
      , PlotOptions[Axes] // Evaluate
    ]
      // Ex["slope_test/wedge_obtuse-slope_test-neighbourhood-slope.png"]
  }
] // AbsoluteTiming


(* ::Subsection:: *)
(*Table for theoretical vs computed corner slopes*)


Module[
  {
    allData,
    apd, alpha,
    gpd, gamma, k,
    theoreticalSlope,
    apdValues, nearCornerFineLengthScaleValues,
    computedSlope, slopeRelError,
    dummyForTrailingCommas
  },
  (* Import all data *)
  allData = Import["slope_test/wedge_obtuse-slope_test-all-data.txt"] // Uncompress;
  allData = Join @@ allData;
  (* Selected value of wedge half-angle *)
  apd = 135;
  alpha = apd * Degree;
  (* Theoretical corner slope *)
  gpd = 60;
  gamma = gpd Degree;
  k = Sin[alpha] / Cos[gamma];
  theoreticalSlope = 1 / Sqrt[k^2 - 1];
  (* Extract computed corner slope data *)
  apdValues = allData[[All, 1]] // Union;
  nearCornerFineLengthScaleValues = allData[[All, 2]] // Union;
  allData // Cases[
    {apd, ncfls_, _, symmetrySlopeData_, __} :> (
      computedSlope[ncfls] = symmetrySlopeData[[1, 2]];
      slopeRelError[ncfls] = computedSlope[ncfls] / theoreticalSlope - 1;
    )
  ];
  (* Build table *)
  Table[
    {ncfls, computedSlope[ncfls], slopeRelError[ncfls]}
    , {ncfls, nearCornerFineLengthScaleValues // Reverse}
  ]
    // TableForm[#,
      , TableHeadings -> {
          None,
          {"ncfls", "computed slope", "rel. error"}
        }
    ] &
] // Ex["wedge_obtuse-slope-comparison-additional-refinement.pdf"]


(* ::Section:: *)
(*Figure: numerical wedge domain (wedge_obtuse-numerical-domain)*)


Module[
  {
    rMax, alpha,
    angleMarkerRadius,
    textStyle,
    dummyForTrailingCommas
  },
  (* Geometry *)
  rMax = rMaxMesh;
  alpha = 135 Degree;
  (* Diagram *)
  angleMarkerRadius = 0.15 rMax;
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  Show[
    (* Domain boundaries *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      (* Wedge walls *)
      Line @ {XYPolar[rMax, -alpha], {0, 0}, XYPolar[rMax, +alpha]},
      (* Far arc at infinity *)
      Circle[{0, 0}, rMax, {-alpha, alpha}],
      {}
    },
    (* Wedge radius *)
    Graphics @ {
      Text[
        rMaxMesh // textStyle
        , XYPolar[rMax / 2, alpha]
        , {2, 0.25}
      ]
    },
    (* Wedge half-angle *)
    Graphics @ {
      Circle[{0, 0}, angleMarkerRadius, {0, alpha}]
    },
    Graphics @ {
      Text[
        "\[Alpha]" // textStyle
        , XYPolar[angleMarkerRadius, alpha/2]
        , {-2.5, -0.6}
      ]
    },
    {}
    , Axes -> True
    , AxesLabel -> Italicise /@ {"x", "y"}
    , ImageSize -> 0.48 ImageSizeTextWidth
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , Method -> {"AxesInFront" -> False}
    , PlotRangePadding -> Scaled[0.067]
    , Ticks -> None
  ]
] // Ex["wedge_obtuse-numerical-domain.pdf"]


(* ::Section:: *)
(*Figure: finite element mesh detail (wedge_obtuse-mesh-detail)*)


Module[
  {
    apd, alpha,
    mesh, meshWireframe,
    rMax, xMin, yMax, xMax,
    tickTextStyle, axisTextStyle,
    dummyForTrailingCommas
  },
  (* Wedge half-angle *)
  apd = 135;
  alpha = apd * Degree;
  (* Import mesh *)
  mesh = Import @ FString["mesh/wedge_obtuse-mesh-apd-{apd}.txt"] // Uncompress // First;
  meshWireframe = mesh["Wireframe"];
  rMax = 0.25;
  {xMin, yMax} = XYPolar[rMax, alpha];
  xMax = Abs[xMin];
  (* Make plot *)
  tickTextStyle = Style[#, LabelSize["Tick"]] & @* LaTeXStyle;
  axisTextStyle = Style[#, LabelSize["Axis"]] & @* LaTeXStyle;
  Show[
    PolarPlot[rMax, {phi, -alpha, alpha}
      , PlotStyle -> None
      , PolarAxes -> {False, True}
      , PolarAxesOrigin -> {alpha, rMax}
      , PolarTicks -> List @
          Table[
            {r, N[r], {0.075 rMax, 0}}
            , {r, FindDivisions[{0, rMax}, 4] // Rest}
          ]
    ]
      /. {Text[Style[r_, {}], coords_, offset_, opts__] :>
        Text[r // tickTextStyle, coords, {2.7, 1.8}, opts]
      }
    ,
    meshWireframe,
    {}
    , ImageSize -> 0.48 * 0.85 ImageSizeTextWidth
    , PlotRange -> {{xMin, xMax}, {-yMax, yMax}}
  ]
] // Ex["wedge_obtuse-mesh-detail.pdf"]


(* ::Section:: *)
(*Figure: numerical wedge solution (wedge_obtuse-solution)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    rMax, xArcCorner, yArcCorner,
    tNumerical,
    x, y, mesh,
    verticalEdge, wallHeight,
    dummyForTrailingCommas
  },
  (* Parameters *)
  apd = 135;
  gpd = 60;
  alpha = apd * Degree;
  gamma = gpd * Degree;
  (* Plot range *)
  rMax = 3;
  {xArcCorner, yArcCorner} = XYPolar[rMax, alpha];
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Plot *)
  mesh = tNumerical["ElementMesh"];
  wallHeight = Ceiling[HHalfPlane[gamma], 0.2];
  verticalEdge @ {x_, y_} :=
    Line @ {{x, y, tNumerical[x, y]}, {x, y, wallHeight}};
  Show[
    (* Numerical solution *)
    Plot3D[
      tNumerical[x, y], Element[{x, y}, mesh]
      , AxesEdge -> {{-1, -1}, {+1, -1}, {-1, -1}}
      , AxesLabel -> {
          Italicise["x"],
          Italicise["y"],
          Italicise["T"] // Margined @ {{0, 5}, {0, 0}}
        }
      , Boxed -> {Back, Bottom, Left}
      , Filling -> 0
      , FillingStyle -> BoundaryTracingStyle["Solution3D"]
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
      , Lighting -> GeneralStyle["AmbientLighting"]
      , PlotPoints -> 20
      , PlotRange -> {0, wallHeight}
      , PlotRangePadding -> {{Scaled[0.125], Automatic}, Automatic, None}
      , PlotStyle -> BoundaryTracingStyle["Solution3D"]
      , TicksStyle -> LabelSize["Tick"]
      , RegionFunction -> Function[{x, y},
          RPolar[x, y] < rMax
        ]
      , ViewPoint -> {3.8, -1.2, 1.7}
    ],
    (* Wedge walls *)
    Plot3D[
      wallHeight, {x, xArcCorner, 0}, {y, -yArcCorner, yArcCorner}
      , Filling -> 0
      , FillingStyle -> BoundaryTracingStyle["Wall3D"]
      , Lighting -> GeneralStyle["AmbientLighting"]
      , Mesh -> None
      , PlotPoints -> 20
      , PlotStyle -> BoundaryTracingStyle["Wall3D"]
      , RegionFunction -> Function[{x, y},
          Abs[y] < x Tan[alpha]
        ]
    ],
    (* Manually drawn solution base edges *)
    Graphics3D @ {
      Line @ {{0, 0, 0}, {xArcCorner, -yArcCorner, 0}},
      {}
    },
    ParametricPlot3D[
      rMax {Cos[phi], Sin[phi], 0}
      , {phi, -alpha, alpha}
      , PlotStyle -> Directive[Black, Thickness[0]]
    ],
    (* Manually drawn wedge wall vertical edges *)
    Graphics3D @ {
      verticalEdge /@ {
        {xArcCorner, -yArcCorner},
        {xArcCorner, +yArcCorner},
        {0, 0}
      },
      {}
    },
    {}
    , ImageSize -> 0.6 ImageSizeTextWidth
  ]
] // Ex["wedge_obtuse-solution.png"
  , Background -> None
  , ImageResolution -> 4 BasicImageResolution
]


(* ::Section:: *)
(*Figure: mesh with additional refinement (wedge_obtuse-mesh-additional-refinement)*)


Module[
  {
    globalCoarseLengthScale, wallFineLengthScale, nearCornerRadius,
    apd, nearCornerFineLengthScale,
      rMax, alpha,
      boundaryPoints, numBoundaryPoints, boundaryElements,
      equilateralTriangleArea, refinementFunction,
      boundaryMesh, mesh, numMeshElements,
      predicateWet,
      tNumerical, tXDerivative, tYDerivative, tSlope,
      rValues,
        symmetryHeightData, symmetrySlopeData,
        wallHeightData, wallSlopeData,
    rMaxPlot, xMinPlot, yMaxPlot, xMaxPlot,
    tickTextStyle, axisTextStyle,
    dummyForTrailingCommas
  },
  (* Constants *)
  globalCoarseLengthScale = 0.2;
  wallFineLengthScale = 0.01;
  nearCornerRadius = 0.01;
  (* Constants (2) *)
  apd = 135;
  nearCornerFineLengthScale = 0.001;
  (* Build mesh *)
    (* Geometry *)
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
  (* Plot *)
  rMaxPlot = 0.055;
  {xMinPlot, yMaxPlot} = XYPolar[rMaxPlot, alpha];
  xMaxPlot = Abs[xMinPlot];
  (* Make plot *)
  tickTextStyle = Style[#, LabelSize["Tick"]] & @* LaTeXStyle;
  axisTextStyle = Style[#, LabelSize["Axis"]] & @* LaTeXStyle;
  Show[
    PolarPlot[rMaxPlot, {phi, -alpha, alpha}
      , PlotStyle -> None
      , PolarAxes -> {False, True}
      , PolarAxesOrigin -> {alpha, rMaxPlot}
      , PolarTicks -> List @
          Table[
            {r, N[r], {0.02, 0}}
            , {r, FindDivisions[{0, rMaxPlot}, 4] // Rest}
          ]
    ]
      /. {Text[Style[r_, {}], coords_, offset_, opts__] :>
        Text[r // tickTextStyle, coords, {2.3, 2.2}, opts]
      }
    ,
    mesh["Wireframe"],
    {}
    , ImageSize -> 0.4 ImageSizeTextWidth
    , PlotRange -> {{xMinPlot, xMaxPlot}, {-yMaxPlot, yMaxPlot}}
  ]
] // Ex["wedge_obtuse-mesh-detail-additional-refinement.pdf"]


(* ::Section:: *)
(*Figure: near-corner slope with additional refinement (wedge_obtuse-slope-additional-refinement-*)*)


Module[
  {
    rMax, casesWithinRange,
    allData,
    slopeMin, slopeMax,
    apd, nearCornerFineLengthScaleValues,
    symmetrySlopeData, wallSlopeData,
    opts,
    symmetryPlot, legend, wallPlot,
    dummyForTrailingCommas
  },
  (* Horizontal plot range *)
  rMax = 0.015;
  casesWithinRange = Cases[{r_, _} /; r <= rMax];
  (* Import all data *)
  allData = Import["slope_test/wedge_obtuse-slope_test-all-data.txt"] // Uncompress;
  allData = Join @@ allData;
  (* Extract data *)
  apd = 135;
  nearCornerFineLengthScaleValues = allData[[All, 2]] // Union;
  allData // Cases[
    {apd, ncfls_, _, ssData_, _, wsData_} :> (
      symmetrySlopeData[ncfls] = ssData // casesWithinRange;
      wallSlopeData[ncfls] = wsData // casesWithinRange;
    )
  ];
  (* Vertical plot range *)
  {slopeMin, slopeMax} = MinMax[
    Table[
      {
        symmetrySlopeData[ncfls][[All, 2]],
        wallSlopeData[ncfls][[All, 2]]
      }
      , {ncfls, nearCornerFineLengthScaleValues}
    ]
  ];
  (* Plot options *)
  opts = {
    AxesLabel -> {Italicise["r"], "Slope"},
    ImageSize -> 0.45 ImageSizeTextWidth,
    Joined -> True,
    LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"],
    TicksStyle -> LabelSize["Tick"],
    PlotMarkers -> {"OpenMarkers", 5},
    PlotRange -> {slopeMin, slopeMax},
    PlotRangeClipping -> False,
    PlotRangePadding -> {{None, 0.001}, {0.05, 0.08}},
    PlotStyle -> Black,
    PlotOptions[Axes] // Evaluate
  };
  (* Slope along line of symmetry *)
  (*
    First use PlotLegends to get a properly styled legend,
    then extract it for separate export.
  *)
  symmetryPlot =
    ListPlot[
      Table[symmetrySlopeData[ncfls], {ncfls, nearCornerFineLengthScaleValues}]
      , PlotLegends -> Placed[
          LineLegend[
            nearCornerFineLengthScaleValues
            , LabelStyle -> LatinModernLabelStyle @ LabelSize["Legend"]
            (* Don't bother with LegendLabel; just write \ell_\ultra in LaTeX *)
            , LegendLayout -> "Row"
            , LegendMarkerSize -> 4 LabelSize["Legend"]
          ]
          , Scaled @ {0.5, 0.5}
        ]
      , opts
    ];
  symmetryPlot /. {
    Legended[p_, Placed[l_, ___]] :> (
      symmetryPlot = p;
      legend = l;
    )
  };
  legend = GraphicsRow[{legend}
    , ImageSize -> ImageSizeTextWidth
    , ItemAspectRatio -> 0.1
    , Spacings -> -ImageSizeTextWidth (* (no idea why this works) *)
  ];
  (* Slope along wall *)
  wallPlot =
    ListPlot[
      Table[wallSlopeData[ncfls], {ncfls, nearCornerFineLengthScaleValues}]
      , opts
    ];
  (* Export *)
  {
    symmetryPlot // Ex["wedge_obtuse-slope-additional-refinement-symmetry.pdf"],
    wallPlot // Ex["wedge_obtuse-slope-additional-refinement-wall.pdf"],
    legend // Ex["wedge_obtuse-slope-additional-refinement-legend.pdf"],
    Nothing
  }
]


(* ::Section:: *)
(*Figure: viable domain  (wedge_obtuse-viable)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical,
    derList, p, q, grad2, f, vi,
    xCritical,
    xMax, xMin, yMax, rMaxWall,
    more, xMinMore, xMaxMore, yMaxMore,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {135, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMax = Ceiling[1.1 xCritical, 0.05];
  xMin = Round[-3 xMax, 0.05];
  yMax = SeekRoot[vi[xMin, #] &, xMin Tan[alpha] {1, 3}, 20];
  rMaxWall = xMin Sec[alpha];
  (* Plot range but more *)
  more = 0.05;
  xMinMore = xMin - more;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Make plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
      , FrameLabel -> {
          Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
          Italicise["y"]
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {
        xMinMore {1, Tan[alpha]},
        {0, 0},
        xMinMore {1, -Tan[alpha]}
      }
    },
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0
      , {x, xMinMore, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 5
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    {}
    , ImageSize -> 0.5 * 0.8 ImageSizeTextWidth
  ]
] // Ex["wedge_obtuse-viable.pdf"]


(* ::Section:: *)
(*Figure: generic traced boundaries  (wedge_obtuse-traced-boundaries)*)


(* (This is a little slow (~5 sec).) *)
Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical,
    derList, p, q, grad2, f, vi,
    xCritical,
    xMax, xMin, yMax, rMaxWall,
    more, xMinMore, xMaxMore, yMaxMore,
    xyStartListWall, xyStartListSymmetry, xyStartList,
    sMax, xyTracedList,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {135, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMax = Ceiling[1.1 xCritical, 0.05];
  xMin = Round[-3 xMax, 0.05];
  yMax = SeekRoot[vi[xMin, #] &, xMin Tan[alpha] {1, 3}, 20];
  rMaxWall = xMin Sec[alpha];
  (* Plot range but more *)
  more = 0.05;
  xMinMore = xMin - more;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Traced boundary starting points (upper branch) *)
  xyStartListWall = Table[XYPolar[r, alpha], {r, Subdivide[rMaxWall, 12]}];
  xyStartListSymmetry = Table[{x, 0}, {x, {0.24, 0.45, 0.65, 0.82} xCritical}];
  xyStartList = Join[xyStartListWall, xyStartListSymmetry];
  (* Traced boundaries (upper branch) *)
  sMax = RPolar[xMax - xMin, 2 yMax];
  xyTracedList =
    Table[
      Quiet[
        ContactTracedBoundary[derList][xyStart, 0, {-1, 1} sMax
          , -1, 1
          , 0
        ]
        , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
      ]
      , {xyStart, xyStartList}
    ];
  (* Make plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
      , FrameLabel -> {
          Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
          Italicise["y"]
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {
        xMinMore {1, Tan[alpha]},
        {0, 0},
        xMinMore {1, -Tan[alpha]}
      }
    },
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0
      , {x, xMinMore, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 5
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Traced boundaries *)
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // IncludeYReflection
          // Evaluate
        , {s, DomainStart[xy], DomainEnd[xy]}
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {xy, xyTracedList}
    ],
    {}
    , ImageSize -> 0.5 * 0.8 ImageSizeTextWidth
  ]
] // Ex["wedge_obtuse-traced-boundaries.pdf"]


(* ::Section:: *)
(*Figure: viable domain, larger scale  (wedge_obtuse-viable-larger-scale)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical,
    derList, p, q, grad2, f, vi,
    xCritical,
    xMax, xMin, yMax, rMaxWall,
    more, xMinMore, xMaxMore, yMaxMore,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {135, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMax = 2 xCritical;
  xMin = -7 xMax;
  yMax = SeekRoot[vi[xMin, #] &, xMin Tan[alpha] {1, 3}, 20];
  rMaxWall = xMin Sec[alpha];
  (* Plot range but more *)
  more = 0.05;
  xMinMore = xMin - more;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Make plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
      , FrameLabel -> {
          Italicise["x"] // Margined @ {{0, 0}, {0, -15}},
          Italicise["y"]
        }
      , FrameTicksStyle -> LabelSize["Tick"]
      , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {
        xMinMore {1, Tan[alpha]},
        {0, 0},
        xMinMore {1, -Tan[alpha]}
      }
    },
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0
      , {x, xMinMore, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 7
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    {}
    , ImageSize -> 0.4 ImageSizeTextWidth
  ]
] // Ex["wedge_obtuse-viable-larger-scale.pdf"]


(* ::Section:: *)
(*Figure: terminal points  (wedge_obtuse-terminal-points)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical,
    derList, p, q, grad2, f, vi,
    xCritical,
    xMax, xMin, yMax, rMaxWall,
    more, xMinMore, xMaxMore, yMaxMore,
    xOrdinary, yOrdinary,
    sMax, xyContourList,
    textStyle, textStyleBracket,
    plot,
    legendLabelStyle,
    legendCurves, legendRegions,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {135, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMax = 3.1 xCritical;
  xMin = -0.65 xMax;
  yMax = SeekRoot[vi[xMin, #] &, xMin Tan[alpha] {1, 3}, 20];
  rMaxWall = xMin Sec[alpha];
  (* Plot range but more *)
  more = 0.05;
  xMinMore = xMin - more;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Ordinary terminal point *)
  xOrdinary = Way[0, xCritical, 0.05];
  yOrdinary = SeekRoot[vi[xOrdinary, #] &, {0, yMax}, 5];
  (* T-contours *)
  sMax = RPolar[xMax - xMin, 2 yMax];
  xyContourList =
    Table[
      ContourByArcLength[tNumerical][xyInitial, 0, sMax {-1, 1}, 1]
      , {xyInitial, {{xCritical, 0}, {xOrdinary, yOrdinary}}}
    ];
  (* Text style *)
  textStyle = Style[#, LabelSize["Point"]] & @* LaTeXStyle;
  textStyleBracket = Style[#, LabelSize["PointBracket"]] &;
  (* Make plot *)
  plot = Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
      , Frame -> None
    ],
    (* Wedge walls *)
    (*Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {
        xMinMore {1, Tan[alpha]},
        {0, 0},
        xMinMore {1, -Tan[alpha]}
      }
    },*)
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0
      , {x, xMinMore, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 5
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Local T-contours *)
    Table[
      ParametricPlot[
        xy[s] // Through
        , {s, DomainStart[xy], DomainEnd[xy]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Background"]
      ]
      , {xy, xyContourList}
    ],
    (* Critical terminal point (x_0, 0) *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ {xCritical, 0}
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
        {xCritical, 0},
        {-1.7, -0.82}
      ] // textStyle,
      Text[
        "critical",
        {xCritical, 0},
        {-1.55, 0.82}
      ] // textStyle,
      {}
    },
    (* Ordinary terminal point *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ {xOrdinary, yOrdinary}
    },
    Graphics @ {
      Text[
        "ordinary",
        {xOrdinary, yOrdinary},
        {-1.4, -0.25}
      ] // textStyle,
      {}
    },
    {}
  ];
  (* Legend *)
  legendLabelStyle = LatinModernLabelStyle @ LabelSize["Legend"];
  legendCurves =
    CurveLegend[
      BoundaryTracingStyle /@ {"Background", "Terminal"},
      {Row @ {Italicise["T"], "\[Hyphen]contour"}, "terminal curve"}
      , LabelStyle -> legendLabelStyle
    ];
  legendRegions =
    RegionLegend[
      BoundaryTracingStyle /@ {"NonViable"},
      {"non\[Hyphen]viable domain"}
      , LabelStyle -> legendLabelStyle
    ];
  (* Combined *)
  GraphicsRow[
    {
      plot,
      Column[Join[legendCurves, legendRegions]
        , Spacings -> {0, {-1.5, -1.5, -1.4}}
      ]
    }
    , ItemAspectRatio -> 2 yMax / (xMax - xMin)
    , ImageSize -> 0.63 ImageSizeTextWidth
    , Spacings -> {0, 0}
  ]
] // Ex["wedge_obtuse-terminal-points.pdf"]


(* ::Section:: *)
(*Figure: generic traced boundary branches  (wedge_obtuse-traced-boundaries-*-branch)*)


(* (This is a little slow (~5 sec).) *)
Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical,
    derList, p, q, grad2, f, vi,
    xCritical,
    xMax, xMin, yMax, rMaxWall,
    more, xMinMore, xMaxMore, yMaxMore,
    rStartListWall, xyStartListWall, xyStartListSymmetry, xyStartList,
    sMax, xyTracedList,
    branchCases, branchYSigns, ySign,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {135, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_obtuse-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMax = 1.2 xCritical;
  xMin = -5 xCritical;
  yMax = SeekRoot[vi[xMin, #] &, xMin Tan[alpha] {1, 3}, 20];
  rMaxWall = xMin Sec[alpha];
  (* Plot range but more *)
  more = 0.05;
  xMinMore = xMin - more;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Traced boundary starting points (upper branch) *)
  rStartListWall = (1/12 + Subdivide[9]) rMaxWall // Prepend[0];
  xyStartListWall = Table[XYPolar[r, alpha], {r, rStartListWall}];
  xyStartListSymmetry = Table[{x, 0}, {x, {0.2, 0.35, 0.5, 0.7, 0.9} xCritical}];
  xyStartList = Join[xyStartListWall, xyStartListSymmetry];
  (* Traced boundaries (upper branch) *)
  sMax = RPolar[xMax - xMin, 2 yMax];
  xyTracedList =
    Table[
      Quiet[
        ContactTracedBoundary[derList][xyStart, 0, {-1, 1} sMax
          , -1, 1
          , 0, Function[{x, y}, x < xMinMore]
        ]
        , {InterpolatingFunction::femdmval, NDSolveValue::nlnum}
      ]
      , {xyStart, xyStartList}
    ];
  (* Make plot *)
  branchCases = {"upper", "lower"};
  branchYSigns = {1, -1};
  Table[
    ySign = AssociationThread[branchCases -> branchYSigns][case];
    Show[
      EmptyFrame[{xMin, xMax}, {-yMax, yMax}
        , Frame -> None
      ],
      (* Wedge walls *)
      Graphics @ {BoundaryTracingStyle["Wall"],
        Line @ {
          xMinMore {1, ySign * Tan[alpha]},
          {0, 0}
        }
      },
      (* Non-viable domain *)
      RegionPlot[
        vi[x, y] < 0
        , {x, xMinMore, xMaxMore}
        , {y, -yMaxMore, yMaxMore}
        , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
        , PlotPoints -> 7
        , PlotStyle -> BoundaryTracingStyle["NonViable"]
      ],
      (* Traced boundaries *)
      Table[
        ParametricPlot[
          EvaluatePair[xy, s, ySign] // Evaluate
          , {s, DomainStart[xy], DomainEnd[xy]}
          , PlotStyle -> BoundaryTracingStyle["Traced"]
        ]
          /. line_Line :> {Arrowheads @ {{-0.08, 0.33}}, Arrow[line]}
        , {xy, xyTracedList}
      ],
      {}
      , ImageSize -> 0.5 * 0.8 ImageSizeTextWidth
    ] // Ex @ FString["wedge_obtuse-traced-boundaries-{case}-branch.pdf"]
    , {case, branchCases}
  ]
]
