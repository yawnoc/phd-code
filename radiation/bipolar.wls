#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*See r3 (manuscripts/radiation-3-bipolar.pdf).*)


(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< Curvilinear`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "bipolar"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Critical terminal points along u = \[Pi]*)


(* ::Subsubsection:: *)
(*(v_\[Natural]\[Pi], A_\[Natural]\[Pi])*)
(*(Warm-to-hot transition hyperbolic coordinate and dimensionless group)*)


(* ::Text:: *)
(*See Section 4.2.1 of r3 (Page r3-8), particularly (r3.43) and (r3.44).*)


vNatPi =
  InverseFunction[
    Function[{v},
      ConditionalExpression[
        Cosh[v] + 1 - v Sinh[v] / 4,
        v > 0
      ]
    ]
  ][0];


aNatPi = Function[{v}, 4 v^3 / Sinh[v]] @ vNatPi;


{vNatPi, aNatPi} // N


(* ::Subsubsection:: *)
(*v_\[Flat]\[Pi] and v_\[Sharp]\[Pi]*)


(* ::Text:: *)
(*See Sections 4.1.1 and 4.1.2 of r3 (Page r3-6).*)


vFlatPi[a_] /; a <= aNatPi :=
  InverseFunction[
    Function[{v},
      ConditionalExpression[
        Cosh[v] + 1 - v^4 / a,
        0 < v <= vNatPi
      ] // Evaluate
    ]
  ][0] // Evaluate;


vSharpPi[a_] /; a <= aNatPi :=
  InverseFunction[
    Function[{v},
      ConditionalExpression[
        Cosh[v] + 1 - v^4 / a,
        v >= vNatPi
      ] // Evaluate
    ]
  ][0] // Evaluate;


(* Check warm-to-hot transition (u = \[Pi]) *)
vFlatPi[aNatPi] == vNatPi == vSharpPi[aNatPi]


(* ::Subsection:: *)
(*Critical terminal points along u = 0*)


(* ::Subsubsection:: *)
(*(v_\[Natural]0, A_\[Natural]0)*)
(*(Cold-to-warm transition hyperbolic coordinate and dimensionless group)*)


(* ::Text:: *)
(*See Section 4.2.2 of r3 (Page r3-8), particularly (r3.49) and (r3.50).*)


vNat0 =
  InverseFunction[
    Function[{v},
      ConditionalExpression[
        Cosh[v] - 1 - v Sinh[v] / 4,
        v > 0
      ]
    ]
  ][0];


aNat0 = Function[{v}, 4 v^3 / Sinh[v]] @ vNat0;


{vNat0, aNat0} // N


(* ::Subsubsection:: *)
(*v_\[Flat]0 and v_\[Sharp]0*)


(* ::Text:: *)
(*See Sections 4.1.1 to 4.1.4 of r3 (Pages r3-6 and r3-7).*)


vFlat0[a_] /; a <= aNat0 :=
  InverseFunction[
    Function[{v},
      ConditionalExpression[
        Cosh[v] - 1 - v^4 / a,
        0 < v <= vNat0
      ] // Evaluate
    ]
  ][0] // Evaluate;


vSharp0[a_] /; a <= aNat0 :=
  InverseFunction[
    Function[{v},
      ConditionalExpression[
        Cosh[v] - 1 - v^4 / a,
        v >= vNat0
      ] // Evaluate
    ]
  ][0] // Evaluate;


(* Check cold-to-warm transition (u = 0) *)
vFlat0[aNat0] == vNat0 == vSharp0[aNat0]


(* ::Subsection:: *)
(*Italicised symbols*)


aIt = Italicised["A"];
uIt = Italicised["u"];
vIt = Italicised["v"];


(* ::Subsection:: *)
(*Options for exported GIFS*)


gifOpts = Sequence[
  AnimationRepetitions -> Infinity,
  "DisplayDurations" -> 0.5
];


(* ::Subsection:: *)
(*Bipolar coordinates*)


(* ::Subsubsection:: *)
(*Effectively v = Infinity*)


vInfinity = 10;


(* ::Subsubsection:: *)
(*Contours*)


defContStyle = LightGray;


uContour[u_, style_: defContStyle] :=
  ParametricPlot[
    XYBipolar[u, v], {v, -vInfinity, vInfinity},
    Axes -> None,
    Exclusions -> If[Mod[u, 2 Pi] == 0, 0, Automatic],
    PlotStyle -> style
  ];


vContour[v_, style_: defContStyle] :=
  ParametricPlot[
    XYBipolar[u, v], {u, 0, 2 Pi},
    Axes -> None,
    PlotStyle -> style
  ];


uContours = uContour /@ Most[Range[0, 360, 30] Degree];
vContours = vContour /@ {-2, -1, -1/2, 0, 1/2, 1, 2};


(* ::Subsubsection:: *)
(*Local orthonormal basis vectors*)


uVector[u_, v_, style_] :=
  Graphics @ {Directive[Thick, style],
    Arrow @ {
      XYBipolar[u, v],
      XYBipolar[u, v] + AUBipolar[u, v]
    }
  };


vVector[u_, v_, style_] :=
  Graphics @ {Directive[Thick, style],
    Arrow @ {
      XYBipolar[u, v],
      XYBipolar[u, v] + AVBipolar[u, v]
    }
  };


(* ::Subsection:: *)
(*Global styles for plots*)


textStyle = Style[#, 18] &;
textStyleBigger = Style[#, 20] &;


natStyle = Red;


flat0Style = Darker[Green];
sharp0Style = Magenta;


flatPiStyle = Blue;
sharpPiStyle = Orange;


nonStyle = Directive[Opacity[0.7], LightGray];
termStyle = Pink;
unphysStyle = Black;


upperStyle = Blue;
lowerStyle = Red;
convexStyle = Black;


guideStyle = Dashed;
pointStyle = PointSize[Large];


uStyle = Blue;
vStyle = Red;


lighter = Directive[Opacity[0.5], #] &;


transStyle =
  Function[{a},
    If[a == aNat0 || a == aNatPi,
      Style[#, natStyle] &,
      Identity
    ]
  ];


(* ::Section:: *)
(*Bipolar coordinates*)


(* ::Subsection:: *)
(*Interactive visualiser*)


DynamicModule[
 {uInit, vInit, xyInit,
  vMax, xMax, yMax
 },
  uInit = 90 Degree;
  vInit = 3/10;
  xyInit = XYBipolar[uInit, vInit];
  vMax = 5;
  xMax = 4;
  yMax = 3;
  Manipulate[
    Show[
      EmptyFrame[{-xMax, xMax}, {-yMax, yMax},
        ImageSize -> 480,
        PlotLabel -> BoxedLabel @ Row[
          {
            Style[uIt == Row @ {u / Degree // N, Degree}, uStyle],
            Style[vIt == N[v], vStyle]
          }, ","
        ]
      ],
      (* All u-contours and v-contours *)
      uContours,
      vContours,
      (* Current u-contour and v-contour *)
      uContour[u, lighter @ uStyle],
      vContour[v, lighter @ vStyle],
      (* Current local orthonormal basis vectors *)
      uVector[u, v, uStyle],
      vVector[u, v, vStyle]
    ]
    (* NOTE: In Version 11, TrackingFunction won't work *)
    (* unless endpoint 2 Pi is numericised *)
    (* https://mathematica.stackexchange.com/a/101695 *)
  , {{u, uInit, uIt}, 0, 2 Pi,
      TrackingFunction -> Function[{uNew},
        u = uNew;
        xy = XYBipolar[u, v];
      ]
  }
  , {{v, vInit, vIt}, -vMax, vMax,
      TrackingFunction -> Function[{vNew},
        v = vNew;
        xy = XYBipolar[u, v];
      ]
  }
  , {{xy, xyInit}, Locator,
      TrackingFunction -> Function[{xyNew},
        xy = xyNew;
        u = UBipolar @@ xy // Mod[#, 2 Pi] &;
        v = VBipolar @@ xy;
      ]
  }]
]


(* ::Subsection:: *)
(*Animation (u varying, v fixed)*)


Module[
 {uValues, v,
  xMax, yMax
 },
  uValues = Subdivide[360, 36] Degree // Most;
  v = 0.3;
  xMax = 4;
  yMax = 3;
  Table[
    Show[
      EmptyFrame[{-xMax, xMax}, {-yMax, yMax},
        ImageSize -> 480,
        PlotLabel -> BoxedLabel @ Row[
          {
            Style[uIt == Row @ {u / Degree // N, Degree}, uStyle],
            Style[vIt == N[v], vStyle]
          }, ","
        ]
      ],
      (* All u-contours and v-contours *)
      uContours,
      vContours,
      (* Current u-contour and v-contour *)
      uContour[u, lighter @ uStyle],
      vContour[v, lighter @ vStyle],
      (* Current local orthonormal basis vectors *)
      uVector[u, v, uStyle],
      vVector[u, v, vStyle]
    ]
  , {u, uValues}]
] // Ex["bipolar-sweep-u.gif", gifOpts]


(* ::Subsection:: *)
(*Animation (u fixed, v varying)*)


Module[
 {u, vValues,
  xMax, yMax
 },
  u = 90 Degree;
  vValues = Subdivide[-4, 4, 20];
  xMax = 4;
  yMax = 3;
  Table[
    Show[
      EmptyFrame[{-xMax, xMax}, {-yMax, yMax},
        ImageSize -> 480,
        PlotLabel -> BoxedLabel @ Row[
          {
            Style[uIt == Row @ {u / Degree // N, Degree}, uStyle],
            Style[vIt == N[v], vStyle]
          }, ","
        ]
      ],
      (* All u-contours and v-contours *)
      uContours,
      vContours,
      (* Current u-contour and v-contour *)
      uContour[u, lighter @ uStyle],
      vContour[v, lighter @ vStyle],
      (* Current local orthonormal basis vectors *)
      uVector[u, v, uStyle],
      vVector[u, v, vStyle]
    ]
  , {v, vValues}]
] // Ex["bipolar-sweep-v.gif", gifOpts]
