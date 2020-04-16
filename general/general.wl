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
(*Figure: Terminal points (terminal-*.pdf)*)


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
    Module[{x, y, s},
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
        (* Local T-contour *)
        ContourPlot[t == 0,
          {x, -xMaxLess, xMaxLess}, {y, -yMaxLess, yMaxLess},
          ContourStyle -> BoundaryTracingStyle["Contour"],
          PlotPoints -> 10
        ],
        (* Non-viable domain and terminal curve *)
        RegionPlot[vi < 0,
          {x, -xMaxMore, xMaxMore}, {y, -yMaxMore, yMaxMore},
          BoundaryStyle -> BoundaryTracingStyle["Terminal"],
          PlotPoints -> 15,
          PlotStyle -> BoundaryTracingStyle["NonViable"]
        ],
        (* Traced boundaries *)
        Table[
          ParametricPlot[xy[s] // Through,
            {s, 1/2 DomainStart[xy], 1/2 DomainEnd[xy]},
            PlotPoints -> 2,
            PlotStyle -> BoundaryTracingStyle["Traced"]
          ]
        , {xy, xyTra}],
        (* Terminal point *)
        Graphics @ {
          GeneralStyle["Point"],
          GeneralStyle["Translucent"],
          Point @ {0, 0}
        }
      ]
    ] // Ex @ FString["terminal-{case}.pdf"]
  , {case, caseList}]
]


(* ::Section:: *)
(*Figure: Curvilinear coordinates (orthogonal-curvilinear-coordinates.pdf)*)


Module[
 {f, fInv, u, v, x, y,
  hu, hv,
  au, av,
  x0, y0, u0, v0,
  du, dv,
  xLeft, xRight, yTop, yBottom,
  contourStyle, originalStyle, displacedStyle,
  orthogonalityMarkerLength, orthogonalityMarkerStyle,
  basisVectorLength, basisVectorStyle,
  textStyle,
  auTipXY, avTipXY,
  huduUV, huduXY,
  hvdvUV, hvdvXY,
  uConstUV, uConstXY,
  vConstUV, vConstXY,
  constText, aText, hText, dText, hdText
 },
  (*
    Any orthogonal curvilinear coordinate system (u, v)
    can be written u + i v == f(x + i y)
    for some analytic function f.
    We take f == exp, f^(-1) == log
    so we have "reversed" polar coordinates.
  *)
  (* Transformation *)
  f = Exp;
  fInv = InverseFunction[f];
  u[x_, y_] := Re @ f[x + I y] // ComplexExpand // Evaluate;
  v[x_, y_] := Im @ f[x + I y] // ComplexExpand // Evaluate;
  x[u_, v_] := Re @ fInv[u + I v] // ComplexExpand // Evaluate;
  y[u_, v_] := (
    Im @ fInv[u + I v]
    // ComplexExpand
    // # /. {Arg[u + I v] -> ArcTan[u, v]} & (* ComplexExpand is dumb *)
    // Evaluate
  );
  (* Local basis *)
  hu[u_, v_] := D[{x, y} @@ {u, v} // Through, u] // Evaluate;
  hv[u_, v_] := D[{x, y} @@ {u, v} // Through, v] // Evaluate;
  (* Local orthonormal basis *)
  au[u_, v_] := Normalize @ hu[u, v];
  av[u_, v_] := Normalize @ hv[u, v];
  (* Original coordinates *)
  {x0, y0} = {0.35, 0.55};
  {u0, v0} = {u, v} @@ {x0, y0} // Through;
  (* Displacements *)
  du = 0.3;
  dv = 0.25;
  (* Plot range *)
  xLeft = 0.25;
  xRight = 0.35;
  yTop = 0.3;
  yBottom = 0.3;
  (* Styles *)
  contourStyle = Black;
  originalStyle = contourStyle;
  displacedStyle = Directive[contourStyle, Dashed];
  orthogonalityMarkerLength = 0.02;
  orthogonalityMarkerStyle = Directive[
    EdgeForm[contourStyle],
    FaceForm[None]
  ];
  basisVectorLength = 0.13;
  basisVectorStyle = Directive[Thickness[0.012], Arrowheads[0.07]];
  textStyle = Style[#, 18] & @* LaTeXStyle;
  (* Coordinates of tips of vectors *)
  auTipXY = {x0, y0} + basisVectorLength au[u0, v0];
  avTipXY = {x0, y0} + basisVectorLength av[u0, v0];
  (* Coordinates of displaced length labels *)
  huduUV = {u0 + du / 2, v0};
  huduXY = {x, y} @@ huduUV // Through;
  hvdvUV = {u0, v0 + dv / 2};
  hvdvXY = {x, y} @@ hvdvUV // Through;
  (* Coordinates of contour labels *)
  uConstUV = {u0, v0 - dv};
  uConstXY = {x, y} @@ uConstUV // Through;
  vConstUV = {u0 - 2/3 du, v0};
  vConstXY = {x, y} @@ vConstUV // Through;
  (* Text functions *)
  constText[coord_] := Italicise[coord] == "const";
  aText[coord_] := Subscript[Embolden["a"], Italicise[coord]];
  hText[coord_] := Subscript[Italicise["h"], Italicise[coord]];
  dText[coord_] :=
    Grid[{{"d", Italicise[coord]}},
      Spacings -> {-0.07, 0}
    ];
  hdText[coord_] :=
    Grid[{{hText[coord], dText[coord]}},
      Spacings -> {0.2, 0}
    ];
  (* Plot *)
  Show[
    (* Original contours *)
    ContourPlot[
      {
        u[x, y] == u0,
        v[x, y] == v0
      },
      {x, x0 - xLeft, x0 + xRight},
      {y, y0 - yBottom, y0 + yTop},
      AspectRatio -> Automatic,
      ContourLabels -> None,
      ContourStyle -> originalStyle,
      Frame -> None,
      ImageSize -> 240
    ],
    Graphics @ {
      (* u == const *)
      Text[
        constText["u"] // textStyle,
        uConstXY,
        {0, -1.2},
        av @@ uConstUV
      ],
      (* v == const *)
      Text[
        constText["v"] // textStyle,
        vConstXY,
        {0, 0.7},
        au @@ vConstUV
      ]
    },
    (* Displaced contours *)
    ContourPlot[
      {
        u[x, y] == u0 + du,
        v[x, y] == v0 + dv
      },
      {x, x0 - xLeft, x0 + xRight},
      {y, y0 - yBottom, y0 + yTop},
      ContourLabels -> None,
      ContourStyle -> displacedStyle
    ],
    (* Orthogonality marker *)
    Graphics @ {orthogonalityMarkerStyle,
      Rotate[
        Rectangle[
          {x0, y0},
          {x0, y0} + orthogonalityMarkerLength
        ],
        ArcTan @@ au[u0, v0],
        {x0, y0}
      ]
    },
    (* Basis vectors *)
    Graphics @ {basisVectorStyle,
      (* a_u *)
      Arrow @ {{x0, y0}, auTipXY},
      Text[
        aText["u"] // textStyle,
        auTipXY,
        {-0.4, -1.6}
      ],
      (* a_v *)
      Arrow @ {{x0, y0}, avTipXY},
      Text[
        aText["v"] // textStyle,
        avTipXY,
        {-1, 1.3}
      ]
    },
    (* Displaced length labels *)
    Graphics @ {
      (* h_u du *)
      Text[
        hdText["u"] // textStyle,
        huduXY,
        {0, 2},
        au @@ huduUV
      ],
      (* h_v du *)
      Text[
        hdText["v"] // textStyle,
        hvdvXY,
        {0, -2.5},
        av @@ hvdvUV
      ]
    }
  ]
] // Ex @ "orthogonal-curvilinear-coordinates.pdf"
