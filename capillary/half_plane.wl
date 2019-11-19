(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< NDSolve`FEM`
<< Conway`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "half_plane"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Finite element mesh*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
Module[{dest},
  dest = "half_plane-mesh.txt";
  ExportIfNotExists[dest,
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
]


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
