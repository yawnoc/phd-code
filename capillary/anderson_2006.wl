(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(*
  This be an independent implementation of
    Anderson M. L., Bassom A. P., & Fowkes N. (2006).
    Exact solutions of the Laplace--Young equation.
    Proc. R. Soc. Lond. Ser. A Math. Phys. Eng. Sci., 462 (2076), 3645--3656.
    <https://doi.org/10.1098/rspa.2006.1744>
  for the purposes of independently reproducing its Figure 3
  (traced boundaries with various indentations for the half-plane solution).
*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< FigureStyles`
SetDirectory @ NotebookDirectory[]


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Known solution*)


(* ::Subsubsection:: *)
(*Half-plane wall height (2.4)*)


etaWall[gamma_] := Sqrt[2 (1 - Sin @ gamma)];


(* ::Subsubsection:: *)
(*Implicit known half-plane solution (2.6)*)


yKnown[eta_] := ArcCosh[2 / eta] - Sqrt[4 - eta^2] - ArcCosh @ Sqrt[2] + Sqrt[2];


(* ::Subsubsection:: *)
(*Solution derivative \[PartialD]\[Eta]/\[PartialD]y*)


etaYDerivative[eta_] := 1 / yKnown'[eta];


(* ::Subsubsubsection:: *)
(*Verify solution derivative satisfies (2.3) with A = 1*)


With[{eta = \[FormalEta], a = 1},
  (1 + etaYDerivative[eta]^2) ^ (-1/2) == a - 1/2 eta^2
    // FullSimplify[#, 0 < eta < Sqrt[2]] &
]


(* ::Subsection:: *)
(*Boundary tracing*)


(* ::Subsubsection:: *)
(*Boundary tracing system of ODEs (3.2)*)


etaWallPlus[gamma_] := Sqrt[2 (1 + Sin @ gamma)];


yDerivative[eta_] :=
  Divide[
    eta^2 - 2,
    eta Sqrt[4 - eta^2]
  ];


(*
  Lower sign is taken for the x-derivative.
  Thus (3.1) should technically have a plus-or-minus sign,
  and (3.2) should have minus-or-plus where it has plus-or-minus.
*)
xDerivative[gamma_][eta_] :=
  Divide[
    -2 yDerivative[eta] Cos[gamma],
    Sqrt[(eta^2 - etaWall[gamma]^2) (etaWallPlus[gamma]^2 - eta^2)]
  ];


(* ::Subsubsubsection:: *)
(*Verify tracing system satisfies boundary condition (3.1)*)


(* Mathematica needs a bit of help *)
assumedPositiveQuantity[gamma_][eta_] := 2 - 4 eta^2 + eta^4 + 2 Cos[2 gamma];


With[{eta = \[FormalEta], gamma = \[FormalGamma]},
  FullSimplify[
    Divide[
      xDerivative[gamma][eta] etaYDerivative[eta],
      Sqrt[xDerivative[gamma][eta]^2 + yDerivative[eta]^2]
    ]
      ==
    Cos[gamma] Sqrt[1 + etaYDerivative[eta]^2]
    , And[
        0 < eta < etaWall[gamma],
        0 < gamma < Pi/2,
        assumedPositiveQuantity[gamma][eta] > 0,
        True
      ]
  ]
]


(* ::Subsubsubsection:: *)
(*Verify assumedPositiveQuantity is actually positive*)


Manipulate[
  Plot[
    assumedPositiveQuantity[gamma][eta]
    , {eta, 0, etaWall[gamma]}
  ]
  , {gamma, 0, Pi/2 - 10^-6}
]


(* ::Subsubsection:: *)
(*Elliptic integral arguments (unnumbered equation after (3.4))*)


curlyPhi[gamma_][eta_] :=
  ArcSin @ Sqrt @ Divide[
    2 + 2 Sin[gamma] - eta^2,
    4 Sin[gamma]
  ];


muPlus[gamma_] := 2 Sin[gamma] / (1 + Sin[gamma]);
muMinus[gamma_] := 2 Sin[gamma] / (1 - Sin[gamma]);


(* ::Subsubsection:: *)
(*Traced boundaries (3.3) & (3.4)*)


(*
  **IMPORTANT**
  1. There is a typo in (3.3). The very last mu_plus should be a mu_minus.
  2. F and Pi are defined as per Mathematica's EllipticF and EllipticPi,
     with the last parameter being m == k^2,
     rather than k used in Gradsteyn & Ryzhik (1980).
*)
xTraced[gamma_][eta_] := (
  Sqrt[2] Cos[gamma] / Sqrt[1 - Sin[gamma]]
  Subtract[
    EllipticF[curlyPhi[gamma][eta], -muMinus[gamma]],
    1 / (1 + Sin[gamma]) EllipticPi[muPlus[gamma], curlyPhi[gamma][eta], -muMinus[gamma]]
  ]
);


(* Same as known solution. *)
yTraced[eta_] := yKnown[eta];


(* ::Subsubsubsection:: *)
(*Verify traced boundaries*)


With[{eta = \[FormalEta], gamma = \[FormalGamma]},
  FullSimplify[
    {
      xTraced[gamma]'[eta] == xDerivative[gamma][eta],
      yTraced'[eta] == yDerivative[eta]
    }
    , And[
        0 < eta < etaWall[gamma],
        0 < gamma < Pi/2,
        assumedPositiveQuantity[gamma][eta] > 0,
        True
      ]
  ]
]
