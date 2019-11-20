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
(*Italicised symbols*)


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
(*Numerical solution*)


Module[{tSolAss, tSol, mesh, gamma},
  tSolAss = Import["half_plane-solution.txt"] // Uncompress;
  Table[
    tSol = tSolAss[gpd];
    mesh = tSol["ElementMesh"];
    gamma = gpd * Degree;
    Block[{x = \[FormalX], y = \[FormalY]},
    (* (Using With results in SetDelayed::wrsym and an empty plot) *)
      Plot3D[tSol[x, y], Element[{x, y}, mesh],
        AxesLabel -> Italicised /@ {"x", "y", "T"},
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
