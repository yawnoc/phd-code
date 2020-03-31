(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< FigureStyles`
SetDirectory @ NotebookDirectory[];


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Section:: *)
(*Terminal points (fig:terminal-*.pdf)*)


(*
  In all cases we put \[CapitalPhi] == y + 3/2 x^2,
  with a terminal point (x, y) == (0, 0).
  The terminal curve is y == -3/2 x^2.
  The non-viable domain is y < -3/2 x^2.
  We choose three known solutions
  with local T-contour T == 0
  to achieve the three (non-degenerate) cases
  for critical terminal points:
    "ordinary"
      T == (x - 1)^2 + (y + 1)^2 - 2
      T == 0 (y ~ x) crosses the terminal curve at a nonzero-angle.
    "critical_hyperbolic"
      T == y + 1/2 x^2
      T == 0 (y ~ 1/2 x^2) is tangential on the viable side.
    "critical_elliptic"
      T == y + 4 x^2
      T == 0 (y ~ -4 x^2) is tangential on the non-viable side.
  Since \[CapitalPhi] == (grad T)^2 - F^2,
  we have F == sqrt((grad T)^2 - \[CapitalPhi]).
 *)


Module[
 {phi,
  caseList, tAss,
  imageSize, xMax, yMax,
  xMaxLess, yMaxLess,
  xMaxMore, yMaxMore,
  vi, t, p, q, grad2, f,
  xyTraSystem,
  sMax, sLow, sHigh, xyTra
 },
  (* Viability function \[CapitalPhi] *)
  phi = Function[{x, y}, y + 3/2 x^2];
  (* List of cases to be plotted *)
  caseList = {
    "ordinary",
    "critical_hyperbolic",
    "critical_elliptic"
  };
  (* Known solutions for the cases *)
  tAss = AssociationThread[
    caseList -> {
      Function[{x, y}, (x - 1)^2 + (y + 1)^2 - 2],
      Function[{x, y}, y + 1/2 x^2],
      Function[{x, y}, y + 4 x^2]
    }
  ];
  (* Plotting constants *)
  imageSize = 180;
  xMax = 0.65;
  yMax = 0.35;
  xMaxLess = 0.95 xMax;
  yMaxLess = 0.95 yMax;
  xMaxMore = 1.2 xMax;
  yMaxMore = 1.2 yMax;
  (* For the three cases *)
  Table[
    With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
      (* Viability function \[CapitalPhi] *)
      vi = phi[x, y];
      (* Known solution T *)
      t = tAss[case][x, y];
      (* Derivatives of T *)
      p = D[t, x];
      q = D[t, y];
      (* Square of gradient *)
      grad2 = p^2 + q^2;
      (* Flux function *)
      f = Sqrt[grad2 - vi];
      (* System of tracing equations *)
      xyTraSystem[upperBranch_] :=
        Module[{sign, xDer, yDer},
          sign = If[upperBranch, +1, -1];
          (* Return system of ODEs *)
          xDer = (-q f + sign p Re @ Sqrt[vi]) / grad2;
          yDer = (+p f + sign q Re @ Sqrt[vi]) / grad2;
          {x' == xDer, y' == yDer } /. {
            x' -> x'[s],
            y' -> y'[s],
            x -> x[s],
            y -> y[s]
          }
        ];
      (* Solve for traced boundaries *)
      If[case != "critical_elliptic",
        (* Not elliptic case: traced boundaries exist *)
        sMax = 1;
        sLow = -sMax;
        sHigh = If[case === "ordinary", 0, sMax];
        xyTra = Table[
          NDSolveValue[
            {
              xyTraSystem[upperBranch],
              x[0] == 0, y[0] == 0
            }, {x, y}, {s, sLow, sHigh},
            NoExtrapolation
          ]
        , {upperBranch, {True, False}}],
        (* Elliptic case: traced boundaries do not exist *)
        xyTra = {}
      ];
      (* Plot *)
      Show[
        EmptyFrame[{-xMax, xMax}, {-yMax, yMax},
          Frame -> None,
          ImageSize -> imageSize
        ],
        (* Non-viable domain and terminal curve *)
        RegionPlot[vi < 0,
          {x, -xMaxMore, xMaxMore}, {y, -yMaxMore, yMaxMore},
          BoundaryStyle -> BoundaryTracingStyle["Terminal"],
          PlotPoints -> 15,
          PlotStyle -> BoundaryTracingStyle["NonViable"]
        ],
        (* Local T-contour *)
        ContourPlot[t == 0,
          {x, -xMaxLess, xMaxLess}, {y, -yMaxLess, yMaxLess},
          ContourStyle -> BoundaryTracingStyle["Contour"],
          PlotPoints -> 10
        ],
        (* Traced boundaries *)
        Table[
          ParametricPlot[xy[s] // Through,
            {s, 1/2 DomainStart[xy], 1/2 DomainEnd[xy]},
            PlotPoints -> 2,
            PlotStyle -> BoundaryTracingStyle["Traced"]
          ]
        , {xy, xyTra}]
      ]
    ] // Ex @ FString["terminal-{case}.pdf"]
  , {case, caseList}]
]
