(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "isoflux"};


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Conformal mapping*)


zOfZeta[zeta_] := zeta^(1/2);
zetaOfZ[z_] := z^2;


(* ::Subsection:: *)
(*Analytic function W*)


wOfZeta[zeta_] := zeta;


(* ::Subsection:: *)
(*Known solution T = Re[W]*)


tOfZeta[zeta_] := Re[wOfZeta[zeta]] // Evaluate;
tOfZ[z_] := tOfZeta[zetaOfZ[z]] // Evaluate;


(* ::Subsection:: *)
(*Slope dW/dz = dW/d\[Zeta] . d\[Zeta]/dz*)


wSlopeOfZ[z_] := wOfZeta'[zetaOfZ[z]] zetaOfZ'[z] // Evaluate;


(* ::Subsection:: *)
(*Flux function f*)


f = 1;


(* ::Subsection:: *)
(*Viability function \[CapitalPhi]*)


phiOfZ[z_] := Abs[wSlopeOfZ[z]]^2 - f // Evaluate;
phiOfZeta[zeta_] := phiOfZ[zOfZeta[zeta]] // Evaluate;


(* ::Subsection:: *)
(*Solution-positive region function*)


solutionPositiveRegionFunction =
  Function[{x, y}, tOfZ[x + I y] > 0 // Evaluate];


(* ::Subsection:: *)
(*Critical terminal point*)


zCritical = Module[{z}, z /. First @ Solve[phiOfZ[z] == 0 && z > 0, z, Reals]];
zetaCritical = zetaOfZ[zCritical];


(* ::Section:: *)
(*Checks*)


zOfZeta[\[FormalZeta]]
zetaOfZ[\[FormalZ]]


tOfZeta[\[FormalZeta]]
tOfZ[\[FormalZ]]


wSlopeOfZ[\[FormalZ]]


f


phiOfZ[\[FormalZ]]
phiOfZeta[\[FormalZeta]]


zCritical
zetaCritical


(* ::Section:: *)
(*Visualisation*)


(* ::Subsection:: *)
(*Known solution T*)


Plot3D[
  tOfZ[x + I y] // Evaluate
  , {x, 0, 3}
  , {y, -3, 3}
  , AspectRatio -> {1, 1, Automatic}
  , PlotRange -> {0, Automatic}
]


(* ::Subsection:: *)
(*Viability function \[CapitalPhi]*)


Plot3D[
  phiOfZ[x + I y] // Evaluate
  , {x, 0, 2}
  , {y, -2, 2}
  , AspectRatio -> {1, 1, Automatic}
  , PlotRange -> {0, Automatic}
  , RegionFunction -> solutionPositiveRegionFunction
]


(* ::Subsection:: *)
(*Non-viable domain & solution contours*)


Show[
  EmptyFrame[
    {0, 1}, {-1, 1}
  ],
  (* Non-viable domain *)
  RegionPlot[
    phiOfZ[x + I y] < 0 && tOfZ[x + I y] > 0
    , {x, 0, 1}
    , {y, -1, 1}
    , BoundaryStyle -> None
    , PlotStyle -> BoundaryTracingStyle["NonViable"]
  ],
  (* Terminal curve *)
  ContourPlot[
    phiOfZ[x + I y] == 0
    , {x, 0, 1}
    , {y, -1, 1}
    , ContourStyle -> BoundaryTracingStyle["Terminal"]
    , RegionFunction -> solutionPositiveRegionFunction
  ],
  (* Solution contours *)
  ContourPlot[
    tOfZ[x + I y]
    , {x, 0, 1}
    , {y, -1, 1}
    , Contours ->
        Module[
          {
            tCriticalValue,
            divisionsUntoCritical,
            divisionsBeyondCritical,
            tSingleDivision,
            dummyForTrailingCommas
          },
          tCriticalValue = tOfZ[zCritical];
          divisionsUntoCritical = 2;
          divisionsBeyondCritical = 5;
          tSingleDivision = tCriticalValue / divisionsUntoCritical;
          tSingleDivision * Range[0, divisionsUntoCritical + divisionsBeyondCritical]
        ]
    , ContourShading -> None
    , ContourLabels -> None
  ],
  {}
]
