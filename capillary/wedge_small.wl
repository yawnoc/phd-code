(* ::Package:: *)

(*
  Basically we extract the singularity
  and treat polar coordinates as "formally rectangular".
  See c1 (manuscripts/capillary-1-small.pdf).
*)


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
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "wedge_small"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Finite element mesh*)


(* Wedge half-angles (apd means "alpha per degree") *)
apdMeshValues = {15, 30, 45, 60};
rMaxMesh = 10;


(* (These are not slow, nevertheless compute once and store.) *)
(* (Delete the files manually to compute from scratch.) *)
Table[
  ExportIfNotExists[FString["mesh/wedge_small-mesh-apd-{apd}.txt"],
    Module[
      {
        rMax, alpha,
        fineLengthScale, coarseLengthScale,
        boundaryPoints, numBoundaryPoints, boundaryElements,
        boundaryMesh, mesh,
        predicateUpperWall, predicateLowerWall, predicateCorner,
        dummyForTrailingCommas
      },
      (* Geometry *)
      rMax = rMaxMesh;
      alpha = apd * Degree;
      (* Refinement length scales *)
      (* (Note that we are using the same refinement for phi as for r) *)
      fineLengthScale = 0.01;
      coarseLengthScale = Min[0.1, alpha/5];
      (* Boundary points *)
      boundaryPoints =
        Join[
          (* Corner segment (r == 0) *)
          Table[{0, phi}, {phi, UniformRange[alpha, -alpha, -fineLengthScale]}] // Most,
          (* Lower wedge wall (\[Phi] == -\[Alpha]) *)
          Table[{r, -alpha}, {r, UniformRange[0, rMaxMesh, fineLengthScale]}] // Most,
          (* Far arc at infinity (r == INF) *)
          Table[{rMaxMesh, phi}, {phi, UniformRange[-alpha, alpha, coarseLengthScale]}] // Most,
          (* Upper wedge wall (\[Phi] == +\[Alpha]) *)
          Table[{r, +alpha}, {r, UniformRange[rMaxMesh, 0, -fineLengthScale]}] // Most,
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
      (* Predicate for wetting boundaries (wedge walls) *)
      predicateUpperWall = Function[{r, phi}, phi == +alpha // Evaluate];
      predicateLowerWall = Function[{r, phi}, phi == -alpha // Evaluate];
      (* Predicate for corner segment (r == 0) *)
      predicateCorner = Function[{r, phi}, r == 0 // Evaluate];
      (* Return *)
      {mesh, predicateUpperWall, predicateLowerWall, predicateCorner} // Compress
    ]
  ]
  , {apd, apdMeshValues}
]


(* ::Subsection:: *)
(*Solve capillary boundary value problem (fixed-point iteration)*)


(* See Page c1-4 of manuscripts/capillary-1-small.pdf *)


(* (This is slow (~3.5 hr), so compute once and store.) *)
(* (Delete the files manually to compute from scratch.) *)
(*
  I do not know how much memory is required for this,
  so I would suggest running this on a fresh kernel.
*)
$HistoryLength = 0;
With[{r = \[FormalR], phi = \[FormalPhi], h = \[FormalCapitalH]},
Module[
  {
    gpdValues,
    mesh, predicateUpperWall, predicateLowerWall, predicateCorner,
    meshCoordinates,
    alpha, gamma, k,
    grad, inactiveGrad, inactiveDiv,
    hCornerSegment, hPrevious,
      nMax, absoluteChangeTolerance,
      pTilde, qTilde, c,
      kMatrix, vVector,
      hCurrent,
      absoluteChange,
      returnList,
    dummyForTrailingCommas
  },
  (* For each value of alpha: *)
  Table[
    (* Solution contact angles (gpd means "gamma per degree") *)
    gpdValues = Join[{1}, Range[5, 90 - apd, 5]];
    (* If all values of contact angle have been account for: *)
    If[
      Table[
        FString @ "solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"
        , {gpd, gpdValues}
      ] // AllTrue[FileExistsQ]
      ,
      (* Then skip; *)
      Null
      ,
      (* Otherwise import mesh *)
      {mesh, predicateUpperWall, predicateLowerWall, predicateCorner} =
        FString["mesh/wedge_small-mesh-apd-{apd}.txt"] // Import // Uncompress;
      meshCoordinates = mesh["Coordinates"];
      (* For each value of gamma: *)
      Table[
        (* Angular parameters etc. *)
        alpha = apd * Degree;
        gamma = gpd * Degree;
        k = Sin[alpha] / Cos[gamma];
        (* Solve and export *)
        ExportIfNotExists[FString @ "solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt",
          (* Formally rectangular derivative operators (see (c1.27)) *)
          grad = Grad[#, {r, phi}] &;
          inactiveGrad = Inactive[Grad][#, {r, phi}] &;
          inactiveDiv = Inactive[Div][#, {r, phi}] &;
          (* Boundary data along corner segment (c1.32) *)
          hCornerSegment = Function[{r, phi}, (Cos[phi] - Sqrt[k^2 - Sin[phi]^2]) / k];
          (* May as well use that as the initial guess *)
          hPrevious = hCornerSegment;
          (* Fixed-point iteration *)
          nMax = 1000;
          absoluteChangeTolerance = 10^-4;
          Do[
            (* Components of gradient de-singularised (see (c1.23) and (c1.24)) *)
            pTilde = r * D[hPrevious[r, phi], r] - hPrevious[r, phi];
            qTilde = D[hPrevious[r, phi], phi];
            (* Common scalar coefficient (c1.26) *)
            c = 1 / Sqrt[r^4 + pTilde^2 + qTilde^2];
            (* Diffusion coefficient matrix (c1.28) *)
            kMatrix = {{r^2 * c, 0}, {0, c}};
            (* Conservative convection coefficients (c1.29) *)
            vVector = {r * c, 0};
            (* Solve linearised boundary value problem in H *)
            hCurrent =
              NDSolveValue[
                {
                  (* Left hand side of (c1.30), negated *)
                  inactiveDiv[
                    -kMatrix . inactiveGrad[h[r, phi]]
                    + Inactive[Times][vVector, h[r, phi]]
                  ] ==
                    (* Right hand side of (c1.30), negated *)
                    - h[r, phi]
                    (* Neumann boundary condition (c1.31) *)
                    + NeumannValue[Cos[gamma], predicateUpperWall[r, phi]]
                    + NeumannValue[Cos[gamma], predicateLowerWall[r, phi]]
                  ,
                  (* Dirichlet condition along corner segment (c1.32) *)
                  DirichletCondition[
                    h[r, phi] == hCornerSegment[r, phi],
                    predicateCorner[r, phi]
                  ]
                },
                h,
                Element[{r, phi}, mesh]
              ];
            (* Compute absolute change *)
            absoluteChange =
              Max @ Table[
                hCurrent @@ coord - hPrevious @@ coord // Abs
                , {coord, meshCoordinates}
              ];
            (* Store return list *)
            returnList = {hCurrent, n, absoluteChange};
            (* If global convergence reached: *)
            If[absoluteChange < absoluteChangeTolerance,
              (* Stop iteration *)
              Break[];
              ,
              (* Otherwise update solution *)
              hPrevious = hCurrent;
            ];
            , {n, nMax}
          ];
          (* Return *)
          returnList // Compress
        ]
        , {gpd, gpdValues}
      ]
    ]
    , {apd, apdMeshValues}
  ]
]
]
$HistoryLength = Infinity;


(* ::Subsection:: *)
(*Boundary tracing*)


(* ::Subsubsection:: *)
(*List of pure-function derivatives*)


(*
  See (c1.23), (c1.24), (c1.37), and (c1.40).
  Analogue to LaplaceYoung`ContactDerivativeList.
*)
tracingDerivativeList[hNumerical_, gammaTracing_?NumericQ] :=
  Module[
    {
      hRDerivative, hPhiDerivative,
      pTilde, qTilde, fTilde, phiTilde,
      dummyForTrailingCommas
    },
    (* Straight derivatives *)
    hRDerivative = Derivative[1, 0][hNumerical];
    hPhiDerivative = Derivative[0, 1][hNumerical];
    (* De-singularised components of gradient (P-tilde, Q-tilde) *)
    pTilde = Function[{r, phi},
      r * hRDerivative[r, phi] - hNumerical[r, phi]
        // Evaluate
    ];
    qTilde = hPhiDerivative;
    (* De-singularised flux function (F-tilde) *)
    fTilde = Function[{r, phi},
      Cos[gammaTracing]
      Sqrt[r^4 + pTilde[r, phi]^2 + qTilde[r, phi]^2]
        // Evaluate
    ];
    (* De-singularised viability function (\[CapitalPhi]-tilde) *)
    phiTilde = Function[{r, phi},
      Sin[gammaTracing]^2 (pTilde[r, phi]^2 + qTilde[r, phi]^2)
        - r^4 * Cos[gammaTracing]^2
        // Evaluate
    ];
    (* Return *)
    {pTilde, qTilde, fTilde, phiTilde}
  ];


(* ::Subsubsection:: *)
(*Elliptic critical terminal point*)


r0[pTilde_, gammaTracing_] :=
  Quiet[
    SeekRoot[
      pTilde[#, 0] / #^2 + Cot[gammaTracing] &,
      {0, rMaxMesh}
    ]
    , {Power::infy}
  ];


(* ::Subsubsection:: *)
(*Traced boundary (arc-length parametrisation)*)


(* See (c1.42) & (c1.43). *)
tracedBoundary[derList_][{r0_, phi0_}, s0_, {sStart_, sEnd_}
  , sSign_: 1
  , branchSign_: 1
  , terminationPhiTilde_: 0
  , terminationFunction_: (False &)
] :=
  Module[
    {
      pTilde, qTilde, fTilde, phiTilde,
      rDerivative, phiDerivative,
      dummyForTrailingCommas
    },
    (* Derivative list *)
    {pTilde, qTilde, fTilde, phiTilde} = derList;
    (* Right hand sides of boundary tracing system of ODEs *)
    rDerivative[r_, phi_] :=
      Divide[
        -qTilde[r, phi] fTilde[r, phi] + branchSign * pTilde[r, phi] Re @ Sqrt @ phiTilde[r, phi],
        pTilde[r, phi]^2 + qTilde[r, phi]^2
      ];
    phiDerivative[r_, phi_] :=
      Divide[
        +pTilde[r, phi] fTilde[r, phi] + branchSign * qTilde[r, phi] Re @ Sqrt @ phiTilde[r, phi],
        pTilde[r, phi]^2 + qTilde[r, phi]^2
      ] / r;
    (* Solve boundary tracing system of ODEs *)
    (* NOTE:
      Using `phi = \[FormalPhi]` results in a Transpose::nmtx error.
      See <https://mathematica.stackexchange.com/q/236730>.
    *)
    With[{r = \[FormalR], phi = \[FormalF], s = \[FormalS]},
      NDSolveValue[
        {
          r'[s] == sSign * rDerivative[r[s], phi[s]],
          phi'[s] == sSign * phiDerivative[r[s], phi[s]],
          r[s0] == r0,
          phi[s0] == phi0,
          WhenEvent[
            {
              phiTilde[r[s], phi[s]] < terminationPhiTilde,
              terminationFunction[r[s], phi[s]]
            },
            "StopIntegration"
          ]
        }
        , {r, phi}
        , {s, sStart, sEnd}
        , NoExtrapolation
      ]
    ]
  ];


(* ::Section:: *)
(*Finite element mesh check*)


Table[
  Module[
    {
      mesh,
        predicateUpperWall, predicateLowerWall, predicateCorner,
      numElements,
      boundaryCoordinates,
        upperWallBoundaryCoordinates,
        lowerWallBoundaryCoordinates,
        cornerBoundaryCoordinates,
      plot,
      dummyForTrailingCommas
    },
    (* Import mesh *)
    {mesh, predicateUpperWall, predicateLowerWall, predicateCorner} =
      FString["mesh/wedge_small-mesh-apd-{apd}.txt"] // Import // Uncompress;
    numElements = Length @ mesh[[2, 1, 1]];
    (* Determine predicate-satisfying boundary coordinates *)
    boundaryCoordinates =
      Part[
        mesh["Coordinates"],
        mesh["BoundaryElements"][[1, 1]] // Flatten // Union
      ];
    upperWallBoundaryCoordinates =
      Select[
        boundaryCoordinates,
        predicateUpperWall @@ # &
      ];
    lowerWallBoundaryCoordinates =
      Select[
        boundaryCoordinates,
        predicateLowerWall @@ # &
      ];
    cornerBoundaryCoordinates =
      Select[
        boundaryCoordinates,
        predicateCorner @@ # &
      ];
    (* Make plot *)
    plot =
      Show[
        mesh["Wireframe"],
        ListPlot[
          {
            upperWallBoundaryCoordinates,
            lowerWallBoundaryCoordinates,
            cornerBoundaryCoordinates
          }
          , PlotMarkers -> {Automatic, Small}
        ],
        {}
        , ImageSize -> Large
        , PlotLabel -> FString @ "{numElements} mesh elements"
      ];
    (* Export *)
    {
      (* Full plot *)
      plot
        // Ex @ FString @ "mesh/wedge_small-mesh-apd-{apd}.png"
      ,
      (* Zoom near origin to verify wetting predicate *)
      Show[plot
        , PlotLabel -> None
        , PlotRange -> {{0, 1}, {1/2, 1}} apd * Degree
      ]
        // Ex @ FString @ "mesh/wedge_small-mesh-apd-{apd}-zoom.png"
    }
  ]
  , {apd, apdMeshValues}
]


(* ::Section:: *)
(*Numerical solution check*)


(* ::Subsection:: *)
(*Table*)


(* (This is slow (~15 sec) from all the calls to Import.) *)
Module[
  {
    gpdValues,
    numIterates, absoluteChange,
    dummyForTrailingCommas
  },
  Join @@ Table[
    gpdValues = Join[{1}, Range[5, 90 - apd, 5]];
    Table[
      {numIterates, absoluteChange} =
        FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
          // Import // Uncompress // Rest;
      {apd * Degree, gpd * Degree, numIterates, absoluteChange}
      , {gpd, gpdValues}
    ]
    , {apd, apdMeshValues}
  ] // TableForm[#
    , TableHeadings -> {None,
        {"alpha", "gamma", "Iterates", "Absolute change"}
      }
  ] &
] // Ex["wedge_small-solution-table.pdf"]


(* ::Subsection:: *)
(*H and T plots*)


(* (This is slow (~30 sec) from all the calls to Import.) *)
Module[
  {
    alpha, gpdValues,
    gamma,
    hSmall, tSmallPolar,
    rMax, phiValues, opts,
    hPlot, tPlot,
    dummyForTrailingCommas
  },
  (* For various alpha *)
  Table[
    alpha = apd * Degree;
    (* For various gamma *)
    gpdValues = Join[{1}, Range[5, 90 - apd, 5]];
    Table[
      gamma = gpd * Degree;
      (* Import wedge_small solution H == H (r, \[Phi]) *)
      hSmall =
        FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
          // Import // Uncompress // First;
      tSmallPolar[r_, phi_] := hSmall[r, phi] / r;
      (* Plotting constants *)
      rMax = rMaxMesh;
      phiValues = Subdivide[0, alpha, 4];
      opts = {
        ImageSize -> 240,
        PlotLabel -> {"apd", apd, "gpd", gpd},
        PlotOptions[Axes]
      };
      (* Make plots *)
      hPlot =
        Plot[
          Table[hSmall[r, phi], {phi, phiValues}] // Evaluate
          , {r, 0, rMax}
          , AxesLabel -> Italicise /@ {"r", "H"}
          , PlotLegends -> LineLegend[phiValues, LegendLabel -> "\[Phi]"]
          , opts // Evaluate
        ];
      tPlot =
        Plot[
          Table[tSmallPolar[r, phi], {phi, phiValues}]
            // Append @ HHalfPlane[gamma]
            // Evaluate
          , {r, 0, rMax}
          , AxesLabel -> Italicise /@ {"r", "T"}
          , PlotStyle -> (
              ConstantArray[Automatic, Length[phiValues]]
                // Append[Dashed]
            )
          , opts // Evaluate
        ];
      (* Export *)
      Row @ {hPlot, tPlot}
        // Ex @ FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.png"]
      , {gpd, gpdValues}
    ]
    , {apd, apdMeshValues}
  ]
]


(* ::Subsection:: *)
(*Comparison with Cartesian solver, borderline case*)


(* (This is slow (~5 sec) from all the calls to Import.) *)
Module[
  {
    alpha, gpdValues,
    gamma,
    hSmall, tSmallPolar,
    tAcute, tAcutePolar,
    rMax, phiValues, opts,
    plot,
    dummyForTrailingCommas
  },
  (* For various alpha *)
  Table[
    alpha = apd * Degree;
    (* For various gamma *)
    gpdValues = {90 - apd};
    Table[
      gamma = gpd * Degree;
      (* Import wedge_small solution H == H (r, \[Phi]) *)
      hSmall =
        FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
          // Import // Uncompress // First;
      tSmallPolar[r_, phi_] := hSmall[r, phi] / r;
      (* Import wedge_acute solution T == T (x, y) *)
      tAcute =
        FString["../wedge_acute/solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
          // Import // Uncompress // First;
      tAcutePolar[r_, phi_] := tAcute @@ XYPolar[r, phi];
      (* Plotting constants *)
      rMax = rMaxMesh;
      phiValues = Subdivide[0, alpha, 4];
      opts = {
        ImageSize -> 270,
        PlotLabel -> {"apd", apd, "gpd", gpd},
        PlotOptions[Axes]
      };
      (* Make plot and export *)
      plot =
        LogPlot[
          Table[
            tSmallPolar[r, phi] / tAcutePolar[r, phi] - 1 // Abs
            , {phi, phiValues}
          ] // Evaluate
          , {r, 0, rMax}
          , AxesLabel -> {Italicise["r"], "Rel. discr."}
          , AxesOrigin -> {0, 0}
          , PlotLegends -> LineLegend[phiValues, LegendLabel -> "\[Phi]"]
          , PlotRange -> Full
          , opts // Evaluate
        ];
      plot
        // Ex @ FString["solution/wedge_small-solution-comparison-apd-{apd}-gpd-{gpd}.png"]
      , {gpd, gpdValues}
    ]
    , {apd, apdMeshValues}
  ]
]


(* ::Section:: *)
(*Figure: original wedge & formal rectangle (wedge_small-domain-*)*)


Module[
  {
    rMax, alpha, numSubdivisions,
    textStyleLabel, opts, paddingLabel, paddingNothing,
    plot,
    xArrow, yArrow,
    dummyForTrailingCommas
  },
  (* Geometry *)
  rMax = 2;
  alpha = 30 Degree;
  numSubdivisions = 4;
  (* Options *)
  textStyleLabel = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  opts = {Nothing
    , ImageSize -> 0.3 ImageSizeTextWidth
    , PlotRange -> {{0, rMax}, rMax Sin[alpha] {-1, 1}}
  };
  paddingLabel = Scaled[0.12];
  paddingNothing = Scaled[0.02];
  (* Original wedge diagram *)
  plot["original-wedge"] = Show[
    Graphics @ {
      (* Constant r (circular arcs) *)
      Table[
        Circle[{0, 0}, r, {-alpha, alpha}]
        , {r, Subdivide[0, rMax, numSubdivisions]}
      ],
      (* Constant \[Phi] (radial rays) *)
      Table[
        Line @ {{0, 0}, XYPolar[rMax, phi]}
        , {phi, Subdivide[-alpha, alpha, numSubdivisions]}
      ],
      (* Labels *)
      Text[
        Italicise["r"] // textStyleLabel
        , XYPolar[rMax / 2, -alpha]
        , {1.7, 0.8}
      ],
      Text[
        "\[Phi]" // textStyleLabel
        , {rMax, 0}
        , {-2.7, 0}
      ],
      {}
    },
    {}
    , PlotRangePadding -> {{paddingNothing, paddingLabel}, paddingNothing}
    , opts
  ];
  (* Formal rectangle diagram *)
  plot["formal-rectangle"] = Show[
    Graphics @ {
      (* Constant r (circular arcs) *)
      Table[
        Line @ {{r, -alpha}, {r, alpha}}
        , {r, Subdivide[0, rMax, numSubdivisions]}
      ],
      (* Constant \[Phi] (radial rays) *)
      Table[
        Line @ {{0, phi}, {rMax, phi}}
        , {phi, Subdivide[-alpha, alpha, numSubdivisions]}
      ],
      (* Labels *)
      Text[
        Italicise["r"] // textStyleLabel
        , {rMax / 2, -alpha}
        , {0, 1}
      ],
      Text[
        "\[Phi]" // textStyleLabel
        , {0, 0}
        , {3, 0}
      ],
      {}
    },
    {}
    , PlotRangePadding -> {{paddingLabel, paddingNothing}, paddingNothing}
    , opts
  ];
  (* Double-arrow between *)
  xArrow = 1;
  yArrow = 0.13;
  plot["double-arrow"] = Show[
    Graphics @ {
      GeneralStyle["VeryThick"], Arrowheads[0.2],
      LightGray,
      Arrow @ {{-xArrow, yArrow}, {xArrow, yArrow}},
      Arrow @ {{xArrow, -yArrow}, {-xArrow, -yArrow}}
    },
    {}
    , ImageSize -> 0.25 ImageSizeTextWidth
    , PlotRange -> All
    , PlotRangePadding -> {yArrow, 2 yArrow}
  ];
  (* Export *)
  Table[
    plot[case]
      // Ex @ FString["wedge_small-domain-{case}.pdf"]
    , {case, {"original-wedge", "formal-rectangle", "double-arrow"}}
  ]
]


(* ::Section:: *)
(*Figure: mesh detail (wedge_small-mesh-detail)*)


Module[
  {
    apd, alpha, rMax,
    mesh, meshWireframe,
    dummyForTrailingCommas
  },
  (* Geometry *)
  apd = 30;
  alpha = apd * Degree;
  rMax = rMaxMesh;
  (* Import mesh *)
  mesh =
    Import @ FString["mesh/wedge_small-mesh-apd-{apd}.txt"]
      // Uncompress // First;
  meshWireframe = mesh["Wireframe"];
  (* Mesh detail *)
  (* (not bothering with full mesh) *)
  Show[
    meshWireframe
    , Frame -> {{True, False}, {True, True}}
    , FrameLabel -> {
        Italicise["r"] // Margined @ {{0, 0}, {0, -0.03 ImageSizeTextWidth}},
        "\[Phi]" // LaTeXStyle
      }
    , FrameTicks -> {Automatic,
        Table[
          {
            phipd * Degree,
            SeparatedRow["VeryThin"][phipd, Style["\[Degree]", Magnification -> 1.2]]
          }
          , {phipd, Subdivide[-apd, apd, 4]}
        ]
      }
    , FrameTicksStyle -> LabelSize["Tick"]
    , ImageSize -> 0.8 ImageSizeTextWidth
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotRange -> {{0, 2}, {-alpha, alpha}}
    , PlotRangeClipping -> True
    , PlotOptions[Frame] // Evaluate
  ]
    // Ex["wedge_small-mesh-detail.pdf"]
]


(* ::Section:: *)
(*Figure: numerical solution (wedge_small-solution-*)*)


Module[
  {
    apd, alpha,
    gpdValues, gpdCases,
    rMax, phiValues,
    opts,
    case,
    hSolution, tSolution,
    rArrow,
    hArrowMin, hArrowMax,
    tArrowMin, tArrowMax,
    arrowStyle, textStyle,
    hPlot, tPlot,
    dummyForTrailingCommas
  },
  (* Wedge half-angle *)
  apd = 30;
  alpha = apd * Degree;
  (* Contact angles *)
  gpdValues = {45, 60};
  gpdCases = {"small", "borderline"};
  (* Plot range *)
  rMax = 6;
  phiValues = Range[0, apd, 10] Degree;
  opts = {Nothing
    , ImagePadding -> {{0.5, 0.5}, {0.5, 0.7}} ImageSizeCentimetre
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotPoints -> 2
    , PlotStyle -> Black
    , TicksStyle -> LabelSize["Tick"]
    , PlotOptions[Axes]
  };
  (* For each solution contact angle: *)
  Table[
    case = AssociationThread[gpdValues -> gpdCases][gpd];
    (* Import numerical solution *)
    hSolution =
      FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
        // Import // Uncompress // First;
    tSolution[r_, phi_] := hSolution[r, phi] / r;
    (* Phi arrow *)
    rArrow = 1.7;
    hArrowMin = hSolution[rArrow, 0] - 0.4;
    hArrowMax = hSolution[rArrow, alpha] + 0.7;
    tArrowMin = tSolution[rArrow, 0] - 0.25;
    tArrowMax = tSolution[rArrow, alpha] + 0.45;
    arrowStyle = Directive[
      GeneralStyle["Translucent"],
      GeneralStyle["Thick"],
      Arrowheads[Medium]
    ];
    textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
    (* Make plots *)
    hPlot =
      Plot[
        Table[hSolution[r, phi], {phi, phiValues}] // Evaluate
        , {r, 0, rMax}
        , AxesLabel -> Italicise /@ {"r", "H"}
        , Epilog -> {
            {arrowStyle,
              Arrow @ {{rArrow, hArrowMin}, {rArrow, hArrowMax}}
            },
            {
              Text["\[Phi]" // textStyle
                , {rArrow, hArrowMax}
                , {0, -1}
              ]
            }
          }
        , PlotRange -> {0, 3}
        , opts // Evaluate
      ];
    tPlot =
      Plot[
        Table[tSolution[r, phi], {phi, phiValues}] // Evaluate
        , {r, 0, rMax}
        , AxesLabel -> Italicise /@ {"r", "T"}
        , Epilog -> {
            {arrowStyle,
              Arrow @ {{rArrow, tArrowMin}, {rArrow, tArrowMax}}
            },
            {
              Text["\[Phi]" // textStyle
                , {rArrow, tArrowMax}
                , {0, -1}
              ]
            }
          }
        , PlotRange -> {0, 2.1}
        , opts // Evaluate
      ];
    GraphicsColumn[{hPlot, tPlot}
      , ImageSize -> 0.45 ImageSizeTextWidth
      , ItemAspectRatio -> 0.75
    ]
      // Ex @ FString["wedge_small-solution-{case}.pdf"]
    , {gpd, gpdValues}
  ]
]


(* ::Section:: *)
(*Figure: borderline comparison (wedge_small-borderline-comparison)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    hSmall, tSmallPolar,
    tAcute, tAcutePolar,
    rMax, rValues,
    relativeDiscrepancyTable,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  apd = 30;
  gpd = 60;
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import wedge_small solution H == H (r, \[Phi]) *)
  hSmall =
    FString["solution/wedge_small-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Import // Uncompress // First;
  tSmallPolar[r_, phi_] := hSmall[r, phi] / r;
  (* Import wedge_acute solution T == T (x, y) *)
  tAcute =
    FString["../wedge_acute/solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Import // Uncompress // First;
  tAcutePolar[r_, phi_] := tAcute @@ XYPolar[r, phi];
  (* Values to be plotted *)
  rMax = 5;
  rValues = Subdivide[0, rMax, 20] // Rest // Prepend[10^-6];
  relativeDiscrepancyTable =
    Table[
      Table[
        {
          r,
          tSmallPolar[r, phi] / tAcutePolar[r, phi] - 1
            // Abs
        }
        , {r, rValues}
      ]
      , {phi, {0, alpha}}
    ];
  (* Make plot and export *)
  ListLogPlot[
    relativeDiscrepancyTable
    , AxesLabel -> {
        Italicise["r"] // Margined @ {{0, 1}, {5, 0}},
        Style["Relative discrepancy", Smaller]
      }
    , ImageSize -> {Automatic, 0.35 ImageSizeTextWidth}
    , Joined -> True
    , LabelStyle -> LatinModernLabelStyle @ LabelSize["Axis"]
    , PlotLegends -> LineLegend[
        {
          LaTeXStyle["\[Phi]"] == 0,
          LaTeXStyle["\[Phi]"] == LaTeXStyle["\[Alpha]"]
        }
        , LabelStyle -> LatinModernLabelStyle @ LabelSize["Legend"]
      ]
    , PlotMarkers -> {{Automatic, 6}, {"OpenMarkers", 5}}
    , PlotStyle -> Black
    , TicksStyle -> LabelSize["Tick"]
    , PlotOptions[Axes] // Evaluate
  ]
    // Ex["wedge_small-borderline-comparison.pdf"]
]
