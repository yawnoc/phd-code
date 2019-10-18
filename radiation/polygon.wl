(* ::Package:: *)

(* ::Text:: *)
(*See r4 (manuscripts/radiation-4-polygon.pdf).*)


(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "polygon"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Schwarz--Christoffel transformation*)


(* ::Subsubsection:: *)
(*h_n(\[Zeta])*)


h[n_][zeta_] := Hypergeometric2F1[1/n, 2/n, 1 + 1/n, zeta ^ n];


(* ::Subsubsection:: *)
(*dz/d\[Zeta]*)


zMapDer[n_][zeta_] := 1 / (h[n][1] (1 - zeta ^ n) ^ (2/n)) // Evaluate;


(* ::Subsubsection:: *)
(*z = z(\[Zeta])*)


zMap[n_][zeta_] := zeta h[n][zeta] / h[n][1] // Evaluate;


(* ::Subsubsection:: *)
(*d\[Zeta]/dz*)


zetaMapDer[n_][zeta_] := h[n][1] (1 - zeta ^ n) ^ (2/n) // Evaluate;


(* ::Subsubsection:: *)
(*(Check these)*)


With[{n = \[FormalN], zeta = \[FormalZeta]},
  {
    zMap[n]'[zeta] == zMapDer[n][zeta] // FullSimplify,
    zMapDer[n][zeta] zetaMapDer[n][zeta] == 1
  }
]


(* ::Subsection:: *)
(*Geometric regions*)


(* ::Subsubsection:: *)
(*Unit disk (\[Zeta]-space)*)


diskPoints[num_] := diskPoints[num] = (
  SeedRandom[0];
  RandomPoint[Disk[], num]
);


(* ::Subsubsection:: *)
(*Regular n-gon (z-space)*)


poly[n_] := RegularPolygon[{1, 2 Pi / n}, n];


polyVertices[n_] := CirclePoints[{1, 2 Pi / n}, n];


polyPoints[n_, num_] := polyPoints[n, num] = (
  SeedRandom[0];
  RandomPoint[poly[n], num]
)


(* ::Subsection:: *)
(*Italicised symbols*)


zIt = Italicised["z"];


(* ::Subsection:: *)
(*Global styles for plots*)


pointStyle = PointSize[Large];


(* ::Section:: *)
(*Schwarz--Christoffel transformation*)


(* ::Subsection:: *)
(*Forward transformation: disk (\[Zeta]-space) to polygon (z-space)*)


(* ::Text:: *)
(*See https://mathematica.stackexchange.com/a/159086*)


Table[
  Module[{rSamp = 4, phiSamp = 4},
    Show[
      ParametricPlot[r Exp[I phi] // zMap[n] // ReIm,
        {r, 0, 1}, {phi, 0, 2 Pi},
        Mesh -> {rSamp - 1, n phiSamp - 1},
        PlotPoints -> {rSamp, n phiSamp + 1},
        PlotOptions[Frame] // Evaluate
      ],
      Graphics @ {Directive[Red, pointStyle],
        Point @ polyVertices[n]
      }
    ]
  ] // Ex @ StringJoin["schwarz-christoffel-forward-", ToString[n], ".pdf"]
, {n, 3, 7}]


(* ::Subsection:: *)
(*Backward transformation: polygon (z-space) to disk (\[Zeta]-space)*)


(* ::Subsubsection:: *)
(*Test method 1: call FindRoot (via SeekRoot)*)


zetaTest1[n_][z0_] := SeekRoot[zMap[n][#] - z0 &, z0];


(* ::Subsubsection:: *)
(*Test method 2: integrate d\[Zeta]/dz from 0 to z*)


zetaTest2[n_][z0_] :=
  With[{zeta = \[FormalW]}, (* NOTE: \[FormalZeta] doesn't work *)
    NDSolveValue[
      {
        (* Integrate d\[Zeta]/dz along straight line from z == 0 to z0 *)
        (* but put z == z_0 * t so that t is real *)
        zeta'[t] == zetaMapDer[n] @ zeta[t] * z0,
        zeta[0] == 0
      }, zeta[1], {t, 0, 1}
    ]
  ];


(* ::Subsubsection:: *)
(*Residuals from disk: \[Zeta]_test(z(\[Zeta])) - \[Zeta]*)


Module[{nValues, zetaValues, zetaTestList, resValues},
  nValues = Range[3, 7];
  zetaValues = Complex @@@ diskPoints[100];
  zetaTestList = {zetaTest1, zetaTest2};
  Table[
    resValues =
      Table[
        zetaTest[n] @ zMap[n] @ zeta - zeta
      , {zeta, zetaValues}] // Abs;
    ListPlot[resValues,
      Joined -> True,
      PlotRange -> All
    ]
  , {n, nValues}, {zetaTest, zetaTestList}]
   // TableForm[#, TableHeadings -> {nValues, zetaTestList}] &
] // Ex["schwarz-christoffel-residuals-disk.pdf"]


(* ::Subsubsection:: *)
(*Residuals from polygon: z(\[Zeta]_test(z)) - z*)


Module[{nValues, zetaTestList, zValues, resValues},
  nValues = Range[3, 7];
  zetaTestList = {zetaTest1, zetaTest2};
  Table[
    zValues = Complex @@@ polyPoints[n, 100];
    resValues =
      Table[
        zMap[n] @ zetaTest[n] @ z - z
      , {z, zValues}] // Abs;
    ListPlot[resValues,
      Joined -> True,
      PlotRange -> All
    ]
  , {n, nValues}, {zetaTest, zetaTestList}]
   // TableForm[#, TableHeadings -> {nValues, zetaTestList}] &
] // Ex["schwarz-christoffel-residuals-poly.pdf"]


(* ::Subsubsection:: *)
(*Conclusion*)


(* ::Text:: *)
(*zetaTest1 has much smaller residuals.*)


(* ::Subsubsection:: *)
(*Timings*)


Module[{nValues, zetaTestList, zValues, timingList},
  nValues = Range[3, 7];
  zetaTestList = {zetaTest1, zetaTest2};
  timingList =
    Transpose @ Table[
      zValues = Complex @@@ polyPoints[n, 1000];
      Table[
        zMap[n] @ zetaTest[n] @ z - z
      , {z, zValues}] // AbsoluteTiming // First
    , {n, nValues}, {zetaTest, zetaTestList}];
  ListPlot[timingList,
    AxesLabel -> {Italicised["n"], "Timing"},
    Joined -> True,
    PlotLegends -> zetaTestList,
    PlotRange -> All,
    PlotOptions[Axes] // Evaluate
  ]
] // Ex["schwarz-christoffel-timings.pdf"]


(* ::Subsubsection:: *)
(*Conclusion*)


(* ::Text:: *)
(*zetaTest1 is faster (and the difference increases with n).*)
