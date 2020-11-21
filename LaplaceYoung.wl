(* ::Package:: *)

(* ::Text:: *)
(*LaplaceYoung.wl (for Wolfram Mathematica)*)


(* ::Text:: *)
(*A Wolfram Language Package (.wl) for numerically solving the Laplace--Young equation.*)
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
  MeshRefinementUniform,
  MeshWireframePositions,
  SolveLaplaceYoung,
  SolveLaplaceYoungFixedPoint
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
(*MeshRefinementUniform*)


MeshRefinementUniform::usage = (
  "MeshRefinementUniform[len]\n"
  <> "Mesh refinement function to ensure that all mesh elements "
  <> "have area not exceeding the equilateral triangle of sidelength len."
);


MeshRefinementUniform[len_] :=
  Function[{vertices, area},
    area > Sqrt[3] / 4 * len^2 // Evaluate
  ];


(* ::Subsubsection:: *)
(*MeshWireframePositions*)


MeshWireframePositions::usage = (
  "MeshWireframePositions[mesh, pattern, containsType (def ContainsAny)]\n"
  <> "Returns list of positions of elements of mesh "
  <> "with at least one coordinate satisfying the pattern pattern. "
  <> "To require all elements satisfying the pattern, "
  <> "use containsType ContainsOnly.\n"
  <> "To be used in mesh visualisation, e.g. mesh[\"Wireframe\"[output]]."
);


MeshWireframePositions[
  mesh_ElementMesh,
  pattern_,
  containsType_: ContainsAny
] :=
  Module[{coordinatesPositions, elementPositions},
    coordinatesPositions = Position[mesh["Coordinates"], pattern] // Flatten;
    elementPositions = (
      Position[
        mesh["MeshElements"][[1, 1]],
        list_List /; containsType[list, coordinatesPositions]
      ]
        // Flatten
        // Map[{1, #} &]
    );
    elementPositions
  ];


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
      (* K == 1 / sqrt(1 + (del T)^2) *)
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
      (* Return solution *)
      tSol
    ]
  ];


(* ::Subsubsection:: *)
(*SolveLaplaceYoungFixedPoint*)


SolveLaplaceYoungFixedPoint::usage = (
  "SolveLaplaceYoungFixedPoint[gamma, mesh, prWet, tol (def 10^-6)]\n"
  <> "Solves the Laplace--Young equation over the finite element mesh mesh "
  <> "with contact angle gamma along the portions of the boundary "
  <> "for which prWet[x, y] is True and tolerance tol.\n"
  <> "Uses a simple fixed-point iteration."
);


SolveLaplaceYoungFixedPoint[
  gamma_?NumericQ,
  mesh_ElementMesh,
  prWet_,
  tol : _?NumericQ : 10^-6
] :=
  With[{x = \[FormalX], y = \[FormalY], t = \[FormalCapitalT]},
    Module[
     {grad, iGrad, iDiv,
      tSol, n, relChange,
      k, tSolNew
     },
      grad = Grad[#, {x, y}] &;
      iGrad = Inactive[Grad][#, {x, y}] &;
      iDiv = Inactive[Div][#, {x, y}] &;
      (* Initialise *)
      tSol = 0 &;
      n = 0;
      relChange = Infinity;
      (* Iterate *)
      Monitor[
        While[Abs[relChange] >= tol,
          (* K == 1 / sqrt(1 + (del T)^2) using previous T *)
          k = 1 / Sqrt[1 + # . #] & @ grad @ tSol[x, y];
          (*
            Solve
              del . [K del T] == T          (interior)
                n . [K del T] == cos(gamma) (boundary)
           *)
          tSolNew = Quiet[
            NDSolveValue[
              {
                iDiv[-k * IdentityMatrix[2] . iGrad @ t[x, y]] ==
                  - t[x, y]
                  + NeumannValue[Cos[gamma], prWet[x, y]]
              }, t, Element[{x, y}, mesh]
            ]
          , {NDSolveValue::femibcnd}];
          n += 1;
          (* Check convergence *)
          relChange = (
            Table[
              Quiet[
                tSolNew @@ xy / tSol @@ xy - 1
              , {Power::infy}]
            , {xy, mesh["Coordinates"]}]
              /. {Indeterminate -> 0}
              // Norm[#, Infinity] &
          );
          (* Update solution *)
          tSol = tSolNew;
        ]
      , {"gamma", gamma, "n", n, "relChange", relChange, "tol", N[tol]}];
      (* Return solution and number of iterates *)
      {tSol, n}
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
