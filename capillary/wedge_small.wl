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
        predicateWet, predicateCorner,
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
      predicateWet = Function[{r, phi}, r > 0 && Abs[phi] == alpha // Evaluate];
      (* Predicate for corner segment (r == 0) *)
      predicateCorner = Function[{r, phi}, r == 0 // Evaluate];
      (* Return *)
      {mesh, predicateWet, predicateCorner} // Compress
    ]
  ]
  , {apd, apdMeshValues}
]


(* ::Section:: *)
(*Finite element mesh check*)


Table[
  Module[
    {
      mesh, predicateWet, predicateCorner, numElements,
      boundaryCoordinates,
      wetBoundaryCoordinates, cornerBoundaryCoordinates,
      plot,
      dummyForTrailingCommas
    },
    (* Import mesh *)
    {mesh, predicateWet, predicateCorner} =
      FString["mesh/wedge_small-mesh-apd-{apd}.txt"] // Import // Uncompress;
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
          {wetBoundaryCoordinates, cornerBoundaryCoordinates}
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
