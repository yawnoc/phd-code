(* ::Package:: *)

(* ::Text:: *)
(*Experimental stuff.*)


(* ::Section:: *)
(*Load packages and set export directory*)


SetDirectory @ NotebookDirectory[];
<< Conway`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "sandbox"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Italicised symbols*)


zIt = Italicised["z"];
wIt = Italicised["w"];


(* ::Section:: *)
(*Schwarz--Christoffel map for regular polygons*)


(* ::Subsection:: *)
(*Theory*)


(* ::Text:: *)
(*The conformal transformation from the unit disk (z-space)*)
(*to a polygon (w-space) with vertices w_1, w_2, ..., w_n*)
(*and interior angles a_1 \[Pi], a_2 \[Pi], ..., a_n \[Pi] is given by*)


(* ::Text:: *)
(*  dw/dz = K (1 - z/z_1)^(a_1-1) (1 - z/z_2)^(a_2-1) ... (1 - z/z_n)^(a_n-1)*)


(* ::Text:: *)
(*for some constant K, where z_1, z_2, ..., z_n are the pre-images*)
(*of w_1, w_2, ..., w_n respectively. Note that*)


(* ::Text:: *)
(*  a_1 \[Pi] + a_2 \[Pi] + ... + a_n \[Pi] = (n - 2) \[Pi],*)


(* ::Text:: *)
(*or*)


(* ::Text:: *)
(*  a_1 + a_2 + ... + a_n = n - 2.*)


(* ::Text:: *)
(*In particular, for a regular n-gon we note the symmetry,*)
(*and take z_1, z_2, ..., z_n to be the n-th roots of unity.*)
(*Since all of the interior angles are equal, we have*)


(* ::Text:: *)
(*  a_1 = a_2 = ... = a_n = 1 - 2/n,*)


(* ::Text:: *)
(*and the transformation becomes*)


(* ::Text:: *)
(*  dw/dz *)
(*  = K (1 - z/z_1)^(-2/n) (1 - z/z_2)^(-2/n) ... (1 - z/z_n)^(-2/n)*)
(*  = K (1 - z^n) ^ (-2/n),*)


(* ::Text:: *)
(*which (choosing w = 0 at z = 0) integrates to give*)


(* ::Text:: *)
(*  w = K z * h_n(z),*)


(* ::Text:: *)
(*where*)


(* ::Text:: *)
(*  h_n(z) = 2F1(1/n, 2/n; 1+1/n; z^n).*)


(* ::Text:: *)
(*Finally, we pick*)


(* ::Text:: *)
(* K = 1 / h_n(1)*)


(* ::Text:: *)
(*so that*)


(* ::Text:: *)
(*  w = z * h_n(z) / h_n(1)*)


(* ::Text:: *)
(*maps z = 1 to w = 1.*)


(* ::Subsection:: *)
(*Forward map*)


(* ::Subsubsection:: *)
(*Definition*)


h[n_][z_] := Hypergeometric2F1[1/n, 2/n, 1 + 1/n, z^n];


wMap[n_][z_] := z h[n][z] / h[n][1] // Evaluate;


(* Check *)
With[{n = \[FormalN], z = \[FormalZ]},
  wMap[n]'[z] == (1 / h[n][1]) (1 - z^n) ^ (-2/n)
   // FullSimplify
]


(* ::Subsubsection:: *)
(*Plot*)


(* ::Text:: *)
(*https://mathematica.stackexchange.com/a/159086*)


Table[
  Module[{rSamp = 4, phiSamp = 4},
    Show[
      ParametricPlot[r Exp[I phi] // wMap[n] // ReIm,
        {r, 0, 1}, {phi, 0, 2 Pi},
        Mesh -> {rSamp - 1, n phiSamp - 1},
        PlotPoints -> {rSamp, n phiSamp + 1}
      ],
      Graphics @ {Directive[Red, PointSize[Large]],
        Point @ CirclePoints[{1, 2 Pi / n}, n]
      },
      PlotOptions[Frame] // Evaluate
    ]
  ] // Ex @ StringJoin["schwarz-christoffel-forward-", ToString[n], ".pdf"]
, {n, 3, 7}]


(* ::Subsection:: *)
(*Inverse map*)


(* ::Subsubsection:: *)
(*Definition*)


(* ::Text:: *)
(*Simply solve w(z) = w for z, with initial guess z = w.*)


zMap[n_][w_] := SeekRoot[wMap[n][#] - w &, w];


(* ::Subsubsection:: *)
(*Check*)


Table[
  Module[{nGon, points, zValues},
    SeedRandom[0];
    nGon = RegularPolygon[{1, 2 Pi / n}, n];
    points = RandomPoint[nGon, 100];
    zValues = Complex @@@ points;
    {
      Show[
        Graphics[nGon],
        ListPlot[points, PlotStyle -> Yellow],
        ImageSize -> 180
      ],
      ListPlot[
        Table[
          zMap[n] @ wMap[n][z] - z
        , {z, zValues}],
        AxesLabel -> None,
        ImageSize -> 270,
        Joined -> True,
        PlotLabel -> Row @ {
          "Residuals ",
          zIt @ wIt[zIt] - zIt
        },
        PlotRange -> All,
        PlotOptions[Axes] // Evaluate
      ]
    } // Row[#, " "] &
  ] // Ex @ StringJoin["schwarz-christoffel-residual-", ToString[n], ".pdf"]
, {n, 3, 7}]
