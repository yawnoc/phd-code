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
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "dip_coating"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Contact angle*)


gamma = 10 Degree;


(* ::Subsection:: *)
(*Domain*)


(* ::Subsubsection:: *)
(*Dipped object (rectangular prism)*)


{prismXMax, prismYMax} = {3, 2} / 2;


(* ::Subsubsection:: *)
(*Vat (cylinder)*)


vatRMax = 5;


(* ::Subsubsection:: *)
(*Mesh for region occupied by liquid*)


(* (This is not slow, nevertheless compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["dip_coating-mesh.txt",
  Module[
    {
      fineLengthScale, coarseLengthScale,
      innerBoundaryPoints, innerBoundaryNumber,
      outerBoundaryPoints, outerBoundaryNumber,
      allBoundaryPoints, allBoundaryElements,
      boundaryMesh, mesh,
      dummyForTrailingCommas
    },
    (* Length scales *)
    fineLengthScale = 0.01;
    coarseLengthScale = 0.5;
    (* Inner boundary (prism) *)
    innerBoundaryPoints =
      Join[
        Table[{prismXMax, y}, {y, UniformRange[-prismYMax, prismYMax, +fineLengthScale]}] // Rest,
        Table[{x, prismYMax}, {x, UniformRange[prismXMax, -prismXMax, -fineLengthScale]}] // Rest,
        Table[{-prismXMax, y}, {y, UniformRange[prismYMax, -prismYMax, -fineLengthScale]}] // Rest,
        Table[{x, -prismYMax}, {x, UniformRange[-prismXMax, prismXMax, +fineLengthScale]}] // Rest,
        {}
      ];
    innerBoundaryNumber = Length[innerBoundaryPoints];
    (* Outer boundary (rim of vat) *)
    outerBoundaryPoints =
      Table[
        XYPolar[vatRMax, phi]
        , {phi, UniformRange[0, 2 Pi, coarseLengthScale / vatRMax]}
      ] // Rest;
    outerBoundaryNumber = Length[outerBoundaryPoints];
    (* Full boundary *)
    allBoundaryPoints = Join[innerBoundaryPoints, outerBoundaryPoints];
    allBoundaryElements =
      LineElement /@ {
        Table[{n, n + 1}, {n, innerBoundaryNumber}]
          // Mod[#, innerBoundaryNumber, 1] &
          ,
        Table[{n, n + 1}, {n, outerBoundaryNumber}]
          // Mod[#, outerBoundaryNumber, 1] &
          // # + innerBoundaryNumber &
          ,
        Nothing
      };
    (* Build mesh *)
    boundaryMesh =
      ToBoundaryMesh[
        "Coordinates" -> allBoundaryPoints,
        "BoundaryElements" -> allBoundaryElements,
        {}
      ];
    mesh =
      ToElementMesh[
        boundaryMesh
        , "ImproveBoundaryPosition" -> True
        , MeshRefinementFunction -> MeshRefinementUniform[coarseLengthScale]
        , "RegionHoles" -> {0, 0}
      ];
    mesh // Compress
  ]
]


(* ::Subsubsection:: *)
(*Predicate for wetting boundary*)


predicateWet =
  Function[{x, y},
    RPolar[x, y] < Way[RPolar[prismXMax, prismYMax], vatRMax] // Evaluate
  ] // Evaluate;


(* ::Subsubsection:: *)
(*Solve PDE*)


(* (This is slow (~5 sec), so compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["dip_coating-solution.txt",
  Module[{mesh, tSol, x, y},
    mesh = Import["dip_coating-mesh.txt"] // Uncompress;
    tSol = SolveLaplaceYoung[10 Degree, mesh, predicateWet];
    tSol // Compress
  ]
]


(* ::Section:: *)
(*Mesh*)


Module[{mesh},
  mesh = Import["dip_coating-mesh.txt"] // Uncompress;
  Show[
    mesh["Wireframe"]
    , PlotLabel -> StringTemplate["`` elements"] @ Length @ mesh[[2, 1, 1]]
  ]
] // Ex["dip_coating-mesh.pdf"]


(* ::Section:: *)
(*Solution*)


Module[{tSol, mesh, x, y},
  tSol = Import["dip_coating-solution.txt"] // Uncompress;
  mesh = tSol["ElementMesh"];
  Plot3D[tSol[x, y], Element[{x, y}, mesh]
    , PlotRange -> Full
  ]
] // Ex["dip_coating-solution.png"]
