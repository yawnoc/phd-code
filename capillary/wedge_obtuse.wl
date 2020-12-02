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


(* ::Subsection:: *)
(*Table for theoretical vs computed corner slopes*)


(* (This is slow (~10 sec) from all the calls to Import.) *)
Module[
  {
    apd, alpha,
    gamma, k,
    infiniteSlopeQ, theoreticalSlope,
    tNumerical, computedSlope, slopeRelError,
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
    computedSlope = -Derivative[1, 0][tNumerical][0, 0];
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
  verticalEdge @ {x_, y_} := Line @ {{x, y, 0}, {x, y, wallHeight}};
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
    (* Manually drawn solution vertical edges *)
    Graphics3D @ {
      Line @ {
        {xArcCorner, -yArcCorner, tNumerical[xArcCorner, -yArcCorner]},
        {xArcCorner, -yArcCorner, 0}
      },
      {}
    },
    (* Manually drawn wedge wall vertical edges *)
    Graphics3D @ {
      verticalEdge /@ {{xArcCorner, -yArcCorner}, {xArcCorner, +yArcCorner}},
      Line @ {{0, 0, tNumerical[0, 0]}, {0, 0, wallHeight}},
      {}
    },
    {}
    , ImageSize -> 0.6 ImageSizeTextWidth
  ]
] // Ex["wedge_obtuse-solution.png"
  , Background -> None
  , ImageResolution -> 4 BasicImageResolution
]
