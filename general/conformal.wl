(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
SetDirectory @ NotebookDirectory[];


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Summary*)


(*
  Here we consider a simple example as proof of concept
  for conformal mappings techniques with boundary tracing:
    n . grad T == -T  (convection)
             T == x   (one-dimensional solution)
  This gives "fins" in the shape of circular arcs:
    y == const \[PlusMinus] sqrt(1 - x^2).
  The steps as given in the thesis appendix are as follows:
    1. Choose a conformal mapping z == z(\[Xi])
    2. Choose an analytic function W == W(\[Xi])
       > We simply choose W(\[Zeta]) == \[Zeta], one-dimensional in \[Zeta]-space.
    3. Write F and \[CapitalPhi] in terms of \[Zeta]
       > We have
       >   F == -T == -Re{\[Zeta]}
       > and
       >   \[CapitalPhi] == |d\[Zeta]/dz|^2 - Re^2{\[Zeta]}.
    4. Compute traced boundaries in \[Zeta]-space
       > We have
       >   d\[Zeta]/ds == -i Re{\[Zeta]} \[PlusMinus] sqrt(|d\[Zeta]/dz|^2 - Re^2{\[Zeta]}).
    5. Map the resulting curves back to z-space
*)


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


(* ::Section:: *)
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
