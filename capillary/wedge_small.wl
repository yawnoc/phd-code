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
