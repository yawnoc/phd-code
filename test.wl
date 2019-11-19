(* ::Package:: *)

(* ::Text:: *)
(*https://mathematica.stackexchange.com/a/155268*)
(*.wl files (or .wls in Version 12) are the way to go for Git.*)
(*I am using .wl for backward compatibility with Version 11,*)
(*and because GitHub currently does not apply syntax colouring for .wls files.*)


SetDirectory @ NotebookDirectory[]


(* ::Section:: *)
(*Testing Conway.wl*)


<< Conway`


Plot[x^2, {x, 0, 1},
  PlotLabel -> BoxedLabel["This is a boxed label"],
  PlotOptions[Axes] // Evaluate
]


DefaultColours


Module[{y, iFun},
  iFun = NDSolveValue[{y'[x] == 0, y[0] == 1}, y, {x, 0, 1}, NoExtrapolation];
  {DomainStart[iFun], DomainEnd[iFun], iFun[1 + DomainEnd[iFun]]}
]


EmptyAxes[{0, 1}, {0, 1}]


EmptyFrame[{0, 1}, {0, 1}]


Module[{var = 0},
  {
    "{{1 + 1}} formats as {1 + 1}.",
    "{{var}} formats as {var}."
  } // FString // TableForm
]


Block[{blockVar = 0},
  "Block variables also work: {blockVar}" // FString
]


Module[{outerVar = 0},
  Module[{innerVar = 1},
    "Nested Modules: outerVar is {outerVar} and innerVar is {innerVar}."
      // FString
  ]
]


"a + beta" // PrettyString["beta" -> "\[Beta]"]


Module[{fun, x},
  fun = Cos[#] - # &;
  x = SeekRoot[fun, {0, 1}];
  {x, fun[x]}
]


Module[{fun, x},
  fun = Cos[#] - # &;
  x = SeekRoot[fun, 0.7];
  {x, fun[x]}
]


SeekRootBisection[Cos[#] - # &, {0, 1}]
SeekRootBisection[Cos[#] - # &, {0, 1}, "ReturnIterations" -> True]


UniformRange[0, 1, 0.099]
UniformRange[0, 2 Pi, 5 Degree]


Way[0, 1]


Way["x", "y", "p"]


Module[{points},
  points = RandomPoint[Circle[], 1000];
  points = points // DeleteNearbyPoints[0.2] // SortByPhi;
  ListPlot[points,
    AspectRatio -> Automatic,
    LabelingFunction -> Function @ Placed[#2[[2]], Center],
    PlotRange -> All,
    PlotStyle -> Yellow
  ]
]


(* ::Section:: *)
(*Testing Curvilinear.wl*)


<< Curvilinear`


(* NOTE: For some reason \[FormalPhi] is not protected in Version 11.0. *)
With[{r = \[FormalR], phi = \[FormalPhi], x = \[FormalX], y = \[FormalY]},
  {
    XPolar[r, phi],
    YPolar[r, phi],
    XYPolar[r, phi],
    RPolar[x, y],
    PhiPolar[x, y],
    RPhiPolar[x, y],
    HRPolar[r, phi],
    HPhiPolar[r, phi],
    ARPolar[r, phi],
    APhiPolar[r, phi],
    APhiPolar[r, phi] == Cross @ ARPolar[r, phi]
  } // TableForm[#, TableDepth -> 1] &
]


With[{u = \[FormalU], v = \[FormalV], x = \[FormalX], y = \[FormalY]},
  {
    XBipolar[u, v],
    YBipolar[u, v],
    XYBipolar[u, v],
    UBipolar[x, y],
    VBipolar[x, y],
    UVBipolar[x, y],
    HBipolar[u, v],
    HUBipolar[u, v],
    HVBipolar[u, v],
    AUBipolar[u, v],
    AVBipolar[u, v],
    AVBipolar[u, v] == Cross @ AUBipolar[u, v] // FullSimplify
  } // TableForm[#, TableDepth -> 1] &
]
