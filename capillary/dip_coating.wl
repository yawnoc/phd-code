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


(* ::Section:: *)
(*Mesh*)


Module[{mesh},
  mesh = Import["dip_coating-mesh.txt"] // Uncompress;
  Show[
    mesh["Wireframe"]
    , PlotLabel -> StringTemplate["`` elements"] @ Length @ mesh[[2, 1, 1]]
  ]
] // Ex["dip_coating-mesh.pdf"]
