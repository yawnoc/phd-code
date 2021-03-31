(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< FigureStyles`
SetDirectory @ NotebookDirectory[];


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Conformal map visualisation*)


conformalVisualise[zOfZeta_] :=
  Module[
    {
      xMin, xMax, yMin, yMax,
      reZetaValues, reZetaMin, reZetaMax,
      imZetaValues, imZetaMin, imZetaMax,
      reZetaValuesImportant, coloursImportant,
      imZetaValuesImportant,
      dummyForTrailingCommas
    },
    (* z-space plot range *)
    {xMin, xMax} = {-2, 2};
    {yMin, yMax} = {-2, 2};
    (* \[Zeta]-space contour values *)
    reZetaValues = Subdivide[-4, 4, 16];
    {reZetaMin, reZetaMax} = MinMax[reZetaValues];
    imZetaValues = Subdivide[-4, 4, 16];
    {imZetaMin, imZetaMax} = MinMax[imZetaValues];
    (* \[Zeta]-space important values *)
    reZetaValuesImportant = {0, 1/2, 1};
    coloursImportant = {Black, Orange, Purple};
    imZetaValuesImportant = {-1, 0, 1};
    (* Make plot *)
    Show[
      EmptyFrame[
        {xMin, xMax}, {yMin, yMax}
        , FrameLabel -> None
        , Background -> White
        , ImageSize -> 240
        , PlotLabel -> zOfZeta[\[FormalZeta]]
      ],
      (* Re{\[Zeta]}-contours *)
      ParametricPlot[
        Table[
          zOfZeta[re + I im] // ReIm
          , {re, reZetaValues}
        ] // Evaluate
        , {im, imZetaMin, imZetaMax}
        , PlotRange -> Full
        , PlotStyle -> LightBlue
      ],
      (* Im{\[Zeta]}-contours *)
      ParametricPlot[
        Table[
          zOfZeta[re + I im] // ReIm
          , {im, imZetaValues}
        ] // Evaluate
        , {re, reZetaMin, reZetaMax}
        , PlotRange -> Full
        , PlotStyle -> LightRed
      ],
      (* Re{\[Zeta]}-contours *)
      ParametricPlot[
        Table[
          zOfZeta[re + I im] // ReIm
          , {re, reZetaValuesImportant}
        ] // Evaluate
        , {im, imZetaMin, imZetaMax}
        , PlotRange -> Full
        , PlotStyle -> (Directive[Thick, #] & /@ coloursImportant)
      ],
      (* Points *)
      ListPlot[
        Table[
          zOfZeta[re + I im] // ReIm
          , {im, imZetaValuesImportant}
          , {re, reZetaValuesImportant}
        ]
        , PlotMarkers -> {"OpenMarkers", Medium}
      ],
      {}
    ]
  ]


(* ::Subsection:: *)
(*Choose cosh (effectively elliptic coordinates)*)


(*
  We consider a simple example as proof of concept
  for conformal mappings techniques with boundary tracing:
    n . grad T == -T^4  (i.e. radiation with A == 1)
  The steps as given in the thesis appendix are as follows:
    1. Choose a conformal mapping z == z(\[Xi])
       > Pick z(\[Zeta]) == cosh(1 - \[Zeta])
    2. Choose an analytic function W == W(\[Xi])
       > We simply choose W(\[Zeta]) == \[Zeta], one-dimensional in \[Zeta]-space.
    3. Write F and \[CapitalPhi] in terms of \[Zeta]
       > We have
       >   F == -T == -Re^4{\[Zeta]}
       > and
       >   \[CapitalPhi] == |csch(1 - \[Zeta])|^2 - Re^8{\[Zeta]}.
    4. Compute traced boundaries in \[Zeta]-space
       > We have
       >   d\[Zeta]/ds == -i Re{\[Zeta]} \[PlusMinus] sqrt(|sech(1 - \[Zeta])|^2 - Re^8{\[Zeta]}).
    5. Map the resulting curves back to z-space
*)


(* ::Subsubsection:: *)
(*Transformation*)


zOfZeta[zeta_] := Cosh[1 - zeta];
zetaOfZ[z_] := 1 - ArcCosh[z];


(* Check *)
zOfZeta @ zetaOfZ[\[FormalZ]] == \[FormalZ]


(* ::Subsubsection:: *)
(*Flux*)


flux[zeta_] := -Re[zeta]^4;


(* ::Subsubsection:: *)
(*Gradient squared*)


gradSquared[zeta_] := 1 / zOfZeta'[zeta] // Abs[#]^2 & // Evaluate;


(* ::Subsubsection:: *)
(*Viability*)


viability[zeta_] := gradSquared[zeta] - flux[zeta]^2 // Evaluate;


(* ::Section::Closed:: *)
(*Conformal map exploration*)


Module[{funList},
  funList = {Nothing
    , Identity
    , -(1 - #)^(1/2) &
    , -(1 - #)^(3/2) &
    , -(1 - #)^2 &
    , Cosh[1 - #] &
  };
  Table[
    conformalVisualise[fun]
    , {fun, funList}
  ]
]


(* ::Section:: *)
(*Choose cosh (effectively elliptic coordinates)*)


(* ::Subsection:: *)
(*Known solution*)


Plot3D[
  Re @ zetaOfZ[x + I y]
  , {x, -2, 2}
  , {y, -2, 2}
  , BoxRatios -> Automatic
  , Exclusions -> None
  , PlotRange -> {0, Full}
]


(* ::Subsection:: *)
(*Viability*)


Plot3D[
  viability @ zetaOfZ[x + I y]
  , {x, -2, 2}
  , {y, -2, 2}
  , Exclusions -> None
  , PlotRange -> {0, Automatic}
]


(* ::Subsection:: *)
(*Viable domain*)


Module[
  {
    xMin, xMax, yMin, yMax,
    xMinMore, xMaxMore, yMinMore, yMaxMore,
    dummyForTrailingCommas
  },
  (* Plot range *)
  {xMin, xMax} = {-2, 2};
  {yMin, yMax} = {-1.5, 1.5};
  (* Plot range but more *)
  {xMinMore, xMaxMore, yMinMore, yMaxMore} =
    1.2 {xMin, xMax, yMin, yMax};
  (* Make plot *)
  Show[
    EmptyFrame[{xMin, xMax}, {yMin, yMax}
    ],
    (* Unphysical region *)
    RegionPlot[
      Re @ zetaOfZ[x + I y] < 0
      , {x, xMinMore, xMaxMore}
      , {y, yMinMore, yMaxMore}
      , BoundaryStyle -> None
      , PlotStyle -> BoundaryTracingStyle["Unphysical"]
    ],
    (* T-contours *)
    ContourPlot[
      Re @ zetaOfZ[x + I y]
      , {x, xMinMore, xMaxMore}
      , {y, yMinMore, yMaxMore}
      , Contours -> Subdivide[0, 1, 4]
      , ContourShading -> None
      , ContourStyle -> Blue
      , Exclusions -> None
      , PlotRange -> {0, 1}
    ],
    (* Non-viable domain (nowhere) *)
    RegionPlot[
      viability @ zetaOfZ[x + I y] < 0
      , {x, xMinMore, xMaxMore}
      , {y, yMinMore, yMaxMore}
      , BoundaryStyle -> None
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    {}
  ]
]
