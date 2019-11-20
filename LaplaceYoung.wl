(* ::Package:: *)

(* ::Text:: *)
(*LaplaceYoung.wl (for Wolfram Mathematica)*)


(* ::Text:: *)
(*A Wolfram Language Package (.wl) for numerically solving the Laplace--Young equation.*)
(*(C) 2019 Conway*)
(*ABSOLUTELY NO WARRANTY, i.e. "GOD SAVE YOU"*)


(* ::Section:: *)
(*Package definitions*)


(* ::Subsection:: *)
(*Start of package*)


BeginPackage["LaplaceYoung`", {"NDSolve`FEM`"}];


(* ::Subsection:: *)
(*Clear existing definitions if any*)


Unprotect["LaplaceYoung`*"];
ClearAll["LaplaceYoung`*"];
ClearAll["LaplaceYoung`*`*"];


(* ::Subsection:: *)
(*Mention non-private symbols*)


{
  HHalfPlane,
  XHalfPlaneUniversal,
  XHalfPlane,
  DHalfPlane,
  SolveLaplaceYoung
};


(* ::Subsection:: *)
(*Private scope*)


(* ::Subsubsection:: *)
(*Start of private scope*)


Begin["`Private`"];


(* ::Subsubsection:: *)
(*HHalfPlane*)


HHalfPlane::usage = (
  "HHalfPlane[gamma]\n"
  <> "Returns the wetting height h == sqrt(2 (1 - sin(gamma))) "
  <> "for the half-plane solution."
);


HHalfPlane[gamma_] := Sqrt[2 (1 - Sin[gamma])];


(* ::Subsubsection:: *)
(*XHalfPlaneUniversal*)


XHalfPlaneUniversal::usage = (
  "XHalfPlaneUniversal[t]\n"
  <> "Returns the half-plane universal curve in implicit form "
  <> "x == arccosh(2 / T) - sqrt(4 - T^2) - arccosh sqrt(2) + sqrt(2)."
);


XHalfPlaneUniversal[t_] :=
  ArcCosh[2 / t] - Sqrt[4 - t^2] - ArcCosh @ Sqrt[2] + Sqrt[2]


(* ::Subsubsection:: *)
(*XHalfPlane*)


XHalfPlane::usage = (
  "XHalfPlane[gamma][t]\n"
  <> "Returns the half-plane solution for contact angle gamma along x == 0 "
  <> "in implicit form x == x (T)."
);


XHalfPlane[gamma_][t_] :=
  XHalfPlaneUniversal[t] - XHalfPlaneUniversal @ HHalfPlane[gamma] // Evaluate;


(* ::Subsubsection:: *)
(*DHalfPlane*)


DHalfPlane::usage = (
  "DHalfPlane[gamma, gammaT]\n"
  <> "Returns the offset distance d(gamma, gammaT) "
  <> "along the half-plane universal curve, between the x == const walls "
  <> "corresponding to contact angles gamma and gammaT."
);


DHalfPlane[gamma_, gammaT_] :=
  XHalfPlaneUniversal @ HHalfPlane[gammaT] - XHalfPlaneUniversal @ HHalfPlane[gamma];


(* ::Subsubsection:: *)
(*SolveLaplaceYoung*)


SolveLaplaceYoung::usage = (
  "SolveLaplaceYoung[gamma, mesh, prWet]\n"
  <> "Solves the Laplace--Young equation over the finite element mesh mesh "
  <> "with contact angle gamma along the portions of the boundary "
  <> "for which prWet[x, y] is True.\n"
  <> "Uses the built-in nonlinear capability in NDSolve`FEM` "
  <> "(Version 12 required)."
);


SolveLaplaceYoung[gamma_?NumericQ, mesh_ElementMesh, prWet_] :=
  With[{x = \[FormalX], y = \[FormalY], t = \[FormalCapitalT]},
    Module[{grad, iGrad, iDiv, k, tSol},
      grad = Grad[#, {x, y}] &;
      iGrad = Inactive[Grad][#, {x, y}] &;
      iDiv = Inactive[Div][#, {x, y}] &;
      (* K == 1 / sqrt(1 + (del u)^2) *)
      k = 1 / Sqrt[1 + # . #] & @ grad @ t[x, y];
      (*
        Solve
          del . [K del T] == T          (interior)
            n . [K del T] == cos(gamma) (boundary)
       *)
      tSol = Quiet[
        NDSolveValue[
          {
            iDiv[-k * IdentityMatrix[2] . iGrad @ t[x, y]] ==
              - t[x, y]
              + NeumannValue[Cos[gamma], prWet[x, y]]
          }, t, Element[{x, y}, mesh]
        ]
      , {NDSolveValue::femibcnd}];
      tSol
    ]
  ];


(* ::Subsubsection:: *)
(*End of private scope*)


End[];


(* ::Subsection:: *)
(*Protect definitions*)


Protect["LaplaceYoung`*"];


(* ::Subsection:: *)
(*End of package*)


EndPackage[];
