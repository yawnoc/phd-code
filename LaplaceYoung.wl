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
  SolveLaplaceYoung
};


(* ::Subsection:: *)
(*Private scope*)


(* ::Subsubsection:: *)
(*Start of private scope*)


Begin["`Private`"];


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
