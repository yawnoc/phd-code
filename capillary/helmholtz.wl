(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< NDSolve`FEM`
<< Conway`
<< Curvilinear`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "helmholtz"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Wedge half-angle*)


(* ::Subsubsection:: *)
(*\[Alpha]*)


alpha = Pi/4;


(* ::Subsection:: *)
(*Known solution*)


(* ::Subsubsection:: *)
(*U*)


uKnown[x_, y_] := Exp[(-x+y)/Sqrt[2]] + Exp[(-x-y)/Sqrt[2]];


(* ::Subsubsection:: *)
(*P = \[PartialD]U/\[PartialD]x, Q = \[PartialD]U/\[PartialD]y*)


p = Derivative[1, 0] @ uKnown;
q = Derivative[0, 1] @ uKnown;


(* ::Subsection:: *)
(*Flux function F*)


f = 1;


(* ::Subsection:: *)
(*Viability \[CapitalPhi]*)


vi[x_, y_] := p[x, y]^2 + q[x, y]^2 - f^2 // Expand // Evaluate


(* ::Section:: *)
(*Algebra*)


(* ::Subsection:: *)
(*Check known solution*)


With[{x = \[FormalX], y = \[FormalY]},
  {
    (* Helmholtz equation *)
    Laplacian[uKnown[x, y], {x, y}] == uKnown[x, y],
    (* Boundary conditions *)
    XYPolar[1, alpha + Pi/2] . Grad[uKnown[x, y], {x, y}] == 1 /. {y -> x Tan[alpha]},
    XYPolar[1, -alpha - Pi/2] . Grad[uKnown[x, y], {x, y}] == 1 /. {y -> -x Tan[alpha]},
    (* Derivatives *)
    p[x, y] == -1/Sqrt[2] (Exp[(-x+y)/Sqrt[2]] + Exp[(-x-y)/Sqrt[2]]),
    q[x, y] ==  1/Sqrt[2] (Exp[(-x+y)/Sqrt[2]] - Exp[(-x-y)/Sqrt[2]]),
    Nothing
  } // FullSimplify
]


(* ::Subsection:: *)
(*Gradient squared*)


With[{x = \[FormalX], y = \[FormalY]},
  p[x, y]^2 + q[x, y]^2 // Expand
]


(* ::Subsection:: *)
(*Check viability*)


With[{x = \[FormalX], y = \[FormalY]},
  vi[x, y] == Exp[-Sqrt[2] (x-y)] + Exp[-Sqrt[2] (x+y)] - 1
    // FullSimplify
]
