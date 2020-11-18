(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< NDSolve`FEM`
<< Conway`
<< Curvilinear`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "helmholtz"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Wedge half-angle*)


(* ::Subsubsection:: *)
(*\[Alpha]*)


alpha = Pi/4;


(* ::Subsection:: *)
(*Known solution*)


(* ::Subsubsection:: *)
(*U*)


uKnown[x_, y_] := Exp[(-x+y)/Sqrt[2]] + Exp[(-x-y)/Sqrt[2]];


(* ::Subsubsection:: *)
(*P = \[PartialD]U/\[PartialD]x, Q = \[PartialD]U/\[PartialD]y*)


p = Derivative[1, 0] @ uKnown;
q = Derivative[0, 1] @ uKnown;


(* ::Subsection:: *)
(*Flux function F*)


f = 1;


(* ::Subsection:: *)
(*Viability \[CapitalPhi]*)


vi[x_, y_] := p[x, y]^2 + q[x, y]^2 - f^2 // Expand // Evaluate


(* ::Subsection:: *)
(*Terminal curve x = x(y)*)


xTerm[y_] := 1/Sqrt[2] * Log[2 Cosh[Sqrt[2] y]];


(* ::Section:: *)
(*Algebra*)


(* ::Subsection:: *)
(*Check known solution*)


With[{x = \[FormalX], y = \[FormalY]},
  {
    (* Helmholtz equation *)
    Laplacian[uKnown[x, y], {x, y}] == uKnown[x, y],
    (* Boundary conditions *)
    XYPolar[1, alpha + Pi/2] . Grad[uKnown[x, y], {x, y}] == 1 /. {y -> x Tan[alpha]},
    XYPolar[1, -alpha - Pi/2] . Grad[uKnown[x, y], {x, y}] == 1 /. {y -> -x Tan[alpha]},
    (* Derivatives *)
    p[x, y] == -1/Sqrt[2] (Exp[(-x+y)/Sqrt[2]] + Exp[(-x-y)/Sqrt[2]]),
    q[x, y] ==  1/Sqrt[2] (Exp[(-x+y)/Sqrt[2]] - Exp[(-x-y)/Sqrt[2]]),
    Nothing
  } // FullSimplify
]


(* ::Subsection:: *)
(*Gradient squared*)


With[{x = \[FormalX], y = \[FormalY]},
  p[x, y]^2 + q[x, y]^2 // Expand
]


(* ::Subsection:: *)
(*Check viability*)


With[{x = \[FormalX], y = \[FormalY]},
  vi[x, y] == Exp[-Sqrt[2] (x-y)] + Exp[-Sqrt[2] (x+y)] - 1
    // FullSimplify
]


(* ::Subsection:: *)
(*Check terminal curve*)


With[{x = \[FormalX], y = \[FormalY]},
  vi[x, y] == 0 /. {x -> xTerm[y]}
    // FullSimplify
]


(* ::Section:: *)
(*Figure: wedge domain (helmholtz-wedge-domain)*)


Module[
  {
    wedgeXMax, wedgeYMax,
    clearanceProportion, xMin, xMax, yMax,
    alpha, angleMarkerRadius,
    normalVectorPosition, normalVectorLength,
    orthogonalityMarkerLength, orthogonalityMarkerStyle,
    wallThicknessCorrection,
    textStyle, testStyleOmega,
    dummyForTrailingCommas
  },
  (* Wedge dimensions *)
  wedgeXMax = 1;
  wedgeYMax = 1;
  (* Plot range *)
  clearanceProportion = 1/10;
  xMin = Way[0, wedgeXMax, -clearanceProportion];
  xMax = Way[0, wedgeXMax, 1 + clearanceProportion];
  yMax = Way[0, wedgeYMax, 1 + clearanceProportion];
  (* Wedge half-angle *)
  alpha = Pi/4;
  angleMarkerRadius = 1/5 wedgeXMax;
  (* Normal vector *)
  normalVectorPosition = 1/2 {wedgeXMax, wedgeYMax};
  normalVectorLength = Sqrt[2]/3 wedgeXMax;
  orthogonalityMarkerLength = 1/18 wedgeXMax;
  orthogonalityMarkerStyle = Directive[EdgeForm[Black], FaceForm[None]];
  wallThicknessCorrection = 0.1;
  (* Diagram *)
  textStyle = Style[#, 24] & @* LaTeXStyle;
  testStyleOmega = Style[#, 30] & @* LaTeXStyle;
  Show[
    EmptyAxes[{xMin, xMax}, {-yMax, yMax}
      , AspectRatio -> Automatic
      , AxesStyle -> Darker[Gray]
      , ImageSize -> 240
      , LabelStyle -> LatinModernLabelStyle[24]
      , Method -> {"AxesInFront" -> False}
      , Ticks -> None
    ],
    (* Wedge half-angle *)
    Graphics @ {
      Circle[{0, 0}, angleMarkerRadius, {0, alpha}]
    },
    Graphics @ {
      Text[
        "\[Alpha]" == SeparatedRow["/"]["\[Pi]", 4] // textStyle
        , XYPolar[angleMarkerRadius, alpha / 2]
        , {-1.2, -0.45}
      ]
    },
    (* Normal vector *)
    Graphics @ {Directive[GeneralStyle["Thick"], Arrowheads[Large]],
      Arrow @ {
        normalVectorPosition,
        normalVectorPosition + XYPolar[normalVectorLength, alpha + Pi/2]
      }
    },
    Graphics @ {
      Text[
        Embolden["n"] // textStyle
        , normalVectorPosition + XYPolar[normalVectorLength, alpha + Pi/2]
        , {0, -1}
      ]
    },
    Graphics @ {orthogonalityMarkerStyle,
      Rectangle[
        normalVectorPosition,
        normalVectorPosition + orthogonalityMarkerLength {1 - wallThicknessCorrection, 1}
      ] // Rotate[#, alpha, normalVectorPosition] &
    },
    (* Wedge domain *)
    Graphics @ {
      Text[
        "\[CapitalOmega]" // testStyleOmega
        , {3/4 wedgeXMax, 1/3 wedgeYMax}
        , {-1.2, -0.5}
      ]
    },
    (* Wedge domain boundary *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {{wedgeXMax, wedgeYMax}, {0, 0}, {wedgeXMax, -wedgeYMax}}
    },
    Graphics @ {
      Text[
        "\[PartialD]\[CapitalOmega]" // testStyleOmega
        , 4/5 {wedgeXMax, wedgeYMax}
        , {1, -1}
      ]
    },
    {}
  ]
] // Ex["helmholtz-wedge-domain.pdf"]
