(* ::Package:: *)

(* ::Text:: *)
(*See r1 (manuscripts/radiation-1-plane.pdf).*)


(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "plane"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Traced boundary y = y(x)*)


(* ::Text:: *)
(*See (r1.17) to (r1.20) (Page r1-2 of manuscripts/radiation-1-plane.pdf).*)


yTraDer = Function[{x}, -x^4 / Sqrt[1 - x^8]];
yTra = Function[{x}, Integrate[yTraDer[x], x] // Evaluate];
yTraMax = Abs @ yTra[1];


With[{x = \[FormalX]},
  {yTraDer[x], yTra[x], yTraMax}
]


(* ::Section:: *)
(*Traced boundary algebra*)


(* ::Subsection:: *)
(*General*)


With[{x = \[FormalX]},
  {
    {"dy/dx", yTraDer[x]},
    {"y", yTra[x]},
    {"y(1)", yTraMax},
    {"y(1)", yTraMax // N}
  } // TableForm
] // Ex["plane-traced-general.pdf"]


(* ::Subsection:: *)
(*Near x = 0*)


With[{x = \[FormalX]},
  Module[{yDerAsy, yAsy},
    yDerAsy = yTraDer[x] + O[x]^12;
    yAsy = Integrate[yDerAsy, x];
    {
      {"dy/dx", yDerAsy},
      {"y", yAsy}
    } // TableForm
  ]
] // Ex["plane-traced-x-near-0.pdf"]


(* ::Subsection:: *)
(*Near x = 1*)


With[{xi = \[FormalXi]},
  Module[{yDerAsy, yAsy},
    yDerAsy = -yTraDer[1 - xi] + O[xi]^(1/2);
    yAsy = Integrate[yDerAsy, xi];
    {
      {"dy/d\[Xi]", yDerAsy},
      {"y", yAsy}
    } // TableForm
  ]
] // Ex["plane-traced-x-near-1.pdf"]


(* ::Section:: *)
(*Traced boundary plot*)


(* ::Subsection:: *)
(*Without known solution*)


Plot[{yTra[x], -yTra[x]}, {x, 0, 1},
  AspectRatio -> Automatic,
  PlotRange -> Full,
  PlotStyle -> Black,
  PlotOptions[Axes] // Evaluate
] // Ex["plane-traced.pdf"]


(* ::Subsection:: *)
(*With known solution*)


Show[
  EmptyAxes[{0, 1}, {-yTraMax, yTraMax}],
  (* Known solution *)
  DensityPlot[x, {x, 0, 1}, {y, -yTraMax, yTraMax},
    ColorFunction -> "TemperatureMap",
    RegionFunction -> Function[{x, y}, Abs[y] < -yTra[x]]
  ],
  (* Traced boundaries *)
  Plot[{yTra[x], -yTra[x]}, {x, 0, 1},
    AspectRatio -> Automatic,
    PlotRange -> Full,
    PlotStyle -> Black,
    PlotOptions[Axes] // Evaluate
  ]
] // Ex["plane-traced-with-known.png"]


(* ::Section:: *)
(*Figure: Traced boundaries, single spike (plane-traced-boundary-spike.pdf)*)


Module[{yMax, yMaxMore, xTerm},
  yMax = Ceiling[yTraMax, 0.1];
  yMaxMore = 1.2 yTraMax;
  xTerm = 1;
  Show[
    (* Pair of traced boundaries forming a spike *)
    Plot[{yTra[x], -yTra[x]}, {x, 0, 1},
      AspectRatio -> Automatic,
      ImageSize -> 240,
      LabelStyle -> LatinModernLabelStyle[12],
      PlotRange -> {-yMax, yMax},
      PlotStyle -> BoundaryTracingStyle["Traced"],
      PlotOptions[Axes] // Evaluate
    ],
    (* Critical terminal curve *)
    Graphics @ {BoundaryTracingStyle["Contour"],
      Line @ {{xTerm, -yMaxMore}, {xTerm, yMaxMore}}
    }
  ]
] // Ex["plane-traced-boundary-spike.pdf"]
