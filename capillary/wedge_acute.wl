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
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "wedge_acute"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Wedge half-angles*)


(* apd stands for "alpha per degree". *)
apdValues = {10, 15, 30, 40, 45, 50, 60, 75};


(* ::Subsection:: *)
(*Contact angles*)


(* gpd stands for "gamma per degree". *)
gpdValues = Join[{1}, Range[5, 85, 5]];


(* ::Subsection:: *)
(*Finite element mesh*)


(* (These are not slow, nevertheless compute once and store.) *)
(* (Delete the files manually to compute from scratch.) *)
Table[
  ExportIfNotExists[FString @ "mesh/wedge_acute-mesh-apd-{apd}.txt",
    Module[
     {rMax, alpha,
      lenCoarse,
      rSpacingWet, rValuesWet,
      phiSpacingFar, phiValuesFar,
      bPointList, nB, mod,
      bMesh, mesh,
      prWet
     },
      (* Geometry *)
      rMax = 10;
      alpha = apd * Degree;
      (* Mesh coarse length scale *)
      lenCoarse = 0.2;
      (* Wet boundary (|\[Phi]| == \[Alpha]) points *)
      (* (fine length scale used here) *)
      rSpacingWet = 0.01;
      rValuesWet = UniformRange[0, rMax, rSpacingWet];
      (* Far boundary (R == R_max) points *)
      phiSpacingFar = lenCoarse / rMax;
      phiValuesFar = UniformRange[-alpha, alpha, phiSpacingFar];
      (* Boundary points *)
      bPointList = Join[
        (* Lower wall: \[Phi] == -\[Alpha], 0 < r < R_max *)
        Table[XYPolar[r, -alpha], {r, rValuesWet}] // Most,
        (* Far arc: r == R_max, -\[Alpha] < \[Phi] < \[Alpha] *)
        Table[XYPolar[rMax, phi], {phi, phiValuesFar}] // Most,
        (* Upper wall: \[Phi] == \[Alpha], 0 < r < R_max *)
        Table[XYPolar[r, alpha], {r, rValuesWet // Reverse}] // Most
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
        "ImproveBoundaryPosition" -> True,
        MeshRefinementFunction -> MeshRefinementUniform[lenCoarse]
      ];
      (* Predicate function for wetting boundaries (|\[Phi]| == \[Alpha]) *)
      prWet = Function[{x, y}, x <= rMax Cos[alpha] // Evaluate];
      {mesh, prWet}
        // Compress
    ]
  ]
, {apd, apdValues}]


(* ::Subsection:: *)
(*Solve PDE (Version 12 required)*)


(* (These are slow (~1 hr), so compute once and store.) *)
(* (Delete the files manually to compute from scratch.) *)
(*
  I do not know how much memory is required for this,
  so I would suggest running this on a fresh kernel.
 *)
(* For each value of alpha: *)
Table[
  Module[{mesh, prWet, gamma, eval},
    (* If all values of gamma have been solved for *)
    If[
      Table[
        FString @ "solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"
      , {gpd, gpdValues}]
        // AllTrue[FileExistsQ],
      (* Do nothing *)
      Null,
      (* Otherwise solve for solutions: *)
      (* Import mesh *)
      {mesh, prWet} =
        Import[
          FString @ "mesh/wedge_acute-mesh-apd-{apd}.txt"
        ] // Uncompress;
      (* For each value of gamma: *)
      Table[
        gamma = gpd * Degree;
        ExportIfNotExists[
          FString @ "solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt",
          eval = EvaluationData @ SolveLaplaceYoung[gamma, mesh, prWet];
          eval /@ {"Result", "Success", "FailureType", "MessagesText", "AbsoluteTiming"}
            // Compress
        ]
      , {gpd, gpdValues}]
    ]
  ]
, {apd, apdValues}]


(* ::Section:: *)
(*Finite element mesh*)


(* For each value of alpha: *)
Table[
  Module[{mesh, prWet, nElem},
    (* Import mesh *)
    {mesh, prWet} =
      Import[
        FString @ "mesh/wedge_acute-mesh-apd-{apd}.txt"
      ] // Uncompress;
    nElem = Length @ mesh[[2, 1, 1]];
    (* Plot *)
    Show[
      mesh["Wireframe"],
      PlotLabel -> FString @ "{nElem} mesh elements",
      PlotOptions[Axes] // Evaluate
    ]
  ] // Ex @ FString @ "mesh/wedge_acute-mesh-apd-{apd}.png"
, {apd, apdValues}]


(* ::Section:: *)
(*Numerical solution*)


(* ::Subsection:: *)
(*Using built-in nonlinear solver*)


(* ::Subsubsection:: *)
(*Table*)


(* (This is slow (~1 min) from all the calls to Import.) *)
Join @@ Table[
  Table[
    Module[{tSol, messages, timing, tCorner},
      (* Import numerical solution, error messages and abs. time *)
      {tSol, messages, timing} = Import[
        FString @ "solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"
      ] // Uncompress // Part[#, {1, 4, 5}] &;
      (* Filter out uninformative messages and make pretty *)
      messages = DeleteCases[
        messages,
        str_ /; StringContainsQ[str, "StringForm"]
      ];
      messages = StringRiffle[messages, "\n"];
      messages = messages /. {"" -> "\[HappySmiley]"};
      (* Compute corner height *)
      tCorner = tSol[0, 0];
      (* Build row of table *)
      {apd * Degree, gpd * Degree, tCorner, timing, messages}
    ]
  , {gpd, gpdValues}]
, {apd, apdValues}] // TableForm[#,
  TableAlignments -> {Left, Center},
  TableDepth -> 2,
  TableHeadings -> {
    None,
    {"\[Alpha]", "\[Gamma]", "T_corner", "Abs. time", "Messages"}
  },
  TableSpacing -> {3, 3}
] & // Ex["wedge_acute-solution-table.pdf"]
