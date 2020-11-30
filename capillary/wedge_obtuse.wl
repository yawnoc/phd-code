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
