(* ::Package:: *)

(* ::Text:: *)
(*https://mathematica.stackexchange.com/a/155268*)
(*.wl files (or .wls in Version 12) are the way to go for Git.*)
(*I am using .wl for backward compatibility with Version 11,*)
(*and because GitHub currently does not apply syntax colouring for .wls files.*)


SetDirectory @ ParentDirectory @ NotebookDirectory[]
<< Conway`
<< Curvilinear`
SetDirectory @ NotebookDirectory[]


(* ::Section:: *)
(*Testing Conway.wl*)


Plot[x^2, {x, 0, 1},
  PlotLabel -> BoxedLabel["This is a boxed label"],
  PlotOptions[Axes] // Evaluate
]


DefaultColours


Module[{y, iFun, iFunReInt},
  iFun = NDSolveValue[
    {y'[x] == 0, y[0] == 1}, y, {x, 0, 1},
    MaxStepFraction -> 1/2000,
    NoExtrapolation
  ];
  iFunReInt = ReInterpolate[iFun];
  {
    DomainStart[iFun],
    DomainEnd[iFun],
    iFun[1 + DomainEnd[iFun]],
    iFun["Grid"] // Length,
    iFunReInt["Grid"] // Length
  }
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


(* Compare displayed and exported LaTeXStyle: *)
Module[{expressions},
  expressions = LaTeXStyle @ Flatten @ {
    $Version,
    Join[
      CharacterRange["0", "9"],
      CharacterRange["A", "Z"],
      CharacterRange["a", "z"]
    ]
      // {Identity, Italicise, Embolden}
      // Through,
    Join[
      CharacterRange["\[CapitalAlpha]", "\[CapitalOmega]"],
      CharacterRange["\[Alpha]", "\[Omega]"]
    ],
    "e" ^ ("i" "\[Pi]") + 1 == 0 // TraditionalForm,
    {"\[PartialD]", "\[Epsilon]", "\[CurlyTheta]", "\[CurlyKappa]", "\[Phi]", "\[CurlyRho]", "\[CurlyPi]", Infinity}
  } // Style[#, 24] &;
  {
    expressions,
    expressions // Ex @ FString[
      "test-latex-{ToLowerCase[$SystemID]}-v{Floor[$VersionNumber]}.pdf"
    ]
  }
]


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


Module[{fun, xMin, xMax, xAny, xFirst},
  fun = Function[{x}, (x - 1/2) Sinc[x]];
  (*
    The roots are 1/2, \[Pi], 2\[Pi], ....
    SeekRootBisection does not always find the first root.
    To that end, use SeekFirstRootBisection.
  *)
  xMin = 0;
  xMax = 5/2 Pi;
  xAny = SeekRootBisection[fun, {xMin, xMax}];
  xFirst = SeekFirstRootBisection[fun, {xMin, xMax}];
  Plot[fun[x], {x, xMin, xMax},
    Epilog -> {
      PointSize[Large],
      Red, Point @ {xAny, 0}, Text["Any", {xAny, 0}, {-1, 2}],
      Blue, Point @ {xFirst, 0}, Text["First", {xFirst, 0}, {-1, 2}]
    }
  ]
]


UniformRange[0, 1, 0.099]
UniformRange[0, 2 Pi, 5 Degree]


Table[
  {x, ToName[x]}
, {x, {0, 0., 1, Pi, 1.*^-10, 3 + 0. I, 3 + 4. I, Infinity}}] // TableForm


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