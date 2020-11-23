(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< NDSolve`FEM`
<< Conway`
<< Curvilinear`
<< LaplaceYoung`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "wedge_acute"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Subsection:: *)
(*Wedge half-angles*)


(* apd stands for "alpha per degree". *)
apdValues = {10, 15, 30, 40, 45, 50, 60, 75};


(* ::Subsection:: *)
(*Contact angles*)


(* gpd stands for "gamma per degree". *)
gpdValues = Join[{1}, Range[5, 85, 5]];


(* ::Subsection:: *)
(*Finite element mesh*)


(* (These are not slow, nevertheless compute once and store.) *)
(* (Delete the files manually to compute from scratch.) *)
Table[
  ExportIfNotExists[FString @ "mesh/wedge_acute-mesh-apd-{apd}.txt",
    Module[
     {rMax, alpha,
      lenCoarse,
      rSpacingWet, rValuesWet,
      phiSpacingFar, phiValuesFar,
      bPointList, nB, mod,
      bMesh, mesh,
      prWet
     },
      (* Geometry *)
      rMax = 10;
      alpha = apd * Degree;
      (* Mesh coarse length scale *)
      lenCoarse = 0.2;
      (* Wet boundary (|\[Phi]| == \[Alpha]) points *)
      (* (fine length scale used here) *)
      rSpacingWet = 0.01;
      rValuesWet = UniformRange[0, rMax, rSpacingWet];
      (* Far boundary (R == R_max) points *)
      phiSpacingFar = lenCoarse / rMax;
      phiValuesFar = UniformRange[-alpha, alpha, phiSpacingFar];
      (* Boundary points *)
      bPointList = Join[
        (* Lower wall: \[Phi] == -\[Alpha], 0 < r < R_max *)
        Table[XYPolar[r, -alpha], {r, rValuesWet}] // Most,
        (* Far arc: r == R_max, -\[Alpha] < \[Phi] < \[Alpha] *)
        Table[XYPolar[rMax, phi], {phi, phiValuesFar}] // Most,
        (* Upper wall: \[Phi] == \[Alpha], 0 < r < R_max *)
        Table[XYPolar[r, alpha], {r, rValuesWet // Reverse}] // Most
      ];
      nB = Length[bPointList];
      mod[n_] := Mod[#, n, 1] &;
      (* Boundary element mesh *)
      bMesh = ToBoundaryMesh[
        "Coordinates" -> bPointList,
        "BoundaryElements" -> {
          LineElement[
            Table[{n, n + 1}, {n, nB}] // mod[nB]
          ]
        }
      ];
      (* Build mesh *)
      mesh = ToElementMesh[bMesh,
        "ImproveBoundaryPosition" -> True,
        MeshRefinementFunction -> MeshRefinementUniform[lenCoarse]
      ];
      (* Predicate function for wetting boundaries (|\[Phi]| == \[Alpha]) *)
      prWet = Function[{x, y}, x <= rMax Cos[alpha] // Evaluate];
      {mesh, prWet}
        // Compress
    ]
  ]
, {apd, apdValues}]


(* ::Subsection:: *)
(*Solve PDE (Version 12 required)*)


(* (These are slow (~1 hr), so compute once and store.) *)
(* (Delete the files manually to compute from scratch.) *)
(*
  I do not know how much memory is required for this,
  so I would suggest running this on a fresh kernel.
 *)
(* For each value of alpha: *)
Table[
  Module[{mesh, prWet, gamma, eval},
    (* If all values of gamma have been solved for *)
    If[
      Table[
        FString @ "solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"
      , {gpd, gpdValues}]
        // AllTrue[FileExistsQ],
      (* Do nothing *)
      Null,
      (* Otherwise solve for solutions: *)
      (* Import mesh *)
      {mesh, prWet} =
        Import[
          FString @ "mesh/wedge_acute-mesh-apd-{apd}.txt"
        ] // Uncompress;
      (* For each value of gamma: *)
      Table[
        gamma = gpd * Degree;
        ExportIfNotExists[
          FString @ "solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt",
          eval = EvaluationData @ SolveLaplaceYoung[gamma, mesh, prWet];
          eval /@ {"Result", "Success", "FailureType", "MessagesText", "AbsoluteTiming"}
            // Compress
        ]
      , {gpd, gpdValues}]
    ]
  ]
, {apd, apdValues}]


(* ::Subsection::Closed:: *)
(*System of ODEs for T-contour (old)*)


(* ::Text:: *)
(*See (57) and (58) (Page 13 of manuscripts/boundary-tracing.pdf)*)
(*for arc length parametrisation of a contour.*)


tContourSystem[tSol_, sSign_: 1] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{t, p, q, slope},
      (* Known solution *)
      t = tSol[x, y];
      (* Components of gradient vector *)
      p = D[t, x];
      q = D[t, y];
      (* Magnitude of gradient vector *)
      slope = Sqrt[p ^ 2 + q ^ 2];
      (* Return system of ODEs *)
      {
        x' == q / slope,
        y' == -p / slope
      } /. {
        x' -> Sign[sSign] x'[s],
        y' -> Sign[sSign] y'[s],
        x -> x[s],
        y -> y[s],
        List -> Sequence
      }
    ]
  ] // Evaluate;


(* ::Subsection::Closed:: *)
(*System of ODEs for terminal curve \[CapitalPhi] = 0 (old)*)


(* ::Text:: *)
(*Note that \[CapitalPhi] = 0 is equivalent to (\[Del]T)^2 = cot(\[Gamma])^2,*)
(*so a \[CapitalPhi]-contour is equivalent to a (\[Del]T)^2-contour.*)


(* vi means \[CapitalPhi]. *)
viContourSystem[tSol_, sSign_: 1] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[
     {t, p, q,
      grad2, grad2P, grad2Q,
      grad2Slope
     },
      (* Known solution *)
      t = tSol[x, y];
      (* Components of gradient vector *)
      p = D[t, x];
      q = D[t, y];
      (* Square of gradient vector *)
      grad2 = p^2 + q^2;
      (* Components of gradient of gradient squared *)
      grad2P = D[grad2, x];
      grad2Q = D[grad2, y];
      (* Magnitude of gradient of gradient squared *)
      grad2Slope = Sqrt[grad2P ^ 2 + grad2Q ^ 2];
      (* Return system of ODEs *)
      {
        x' == grad2Q / grad2Slope,
        y' == -grad2P / grad2Slope
      } /. {
        x' -> Sign[sSign] x'[s],
        y' -> Sign[sSign] y'[s],
        x -> x[s],
        y -> y[s],
        List -> Sequence
      }
    ]
  ] // Evaluate;


(* ::Subsection::Closed:: *)
(*System of ODEs for tracing (old)*)


(* ::Text:: *)
(*See (55) and (56) (Page 13 of manuscripts/boundary-tracing.pdf)*)
(*for arc length parametrisation of a traced boundary.*)


xyTraSystem[tSol_, gammaTra_, sSign_: 1] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{t, p, q, grad2, f, vi},
      (* Known solution *)
      t = tSol[x, y];
      (* Components of gradient vector *)
      p = D[t, x];
      q = D[t, y];
      (* Square of gradient vector *)
      grad2 = p^2 + q^2;
      (* Flux function F *)
      f = Cos[gammaTra] Sqrt[1 + grad2];
      (* Viability function \[CapitalPhi] *)
      vi = Sin[gammaTra]^2 * grad2 - Cos[gammaTra]^2;
      (* Return system of ODEs *)
      {
        x' == (-q f + p * Re @ Sqrt[vi]) / grad2,
        y' == (+p f + q * Re @ Sqrt[vi]) / grad2
      } /. {
        x' -> Sign[sSign] x'[s],
        y' -> Sign[sSign] y'[s],
        x -> x[s],
        y -> y[s],
        List -> Sequence
      }
    ]
  ] // Evaluate;


(* ::Text:: *)
(*Termination at terminal curve*)


xyTraSystemTerm[tSol_, gammaTra_] :=
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{t, p, q, grad2, vi},
      (* Known solution *)
      t = tSol[x, y];
      (* Components of gradient vector *)
      p = D[t, x];
      q = D[t, y];
      (* Square of gradient vector *)
      grad2 = p^2 + q^2;
      (* Viability function \[CapitalPhi] *)
      vi = Sin[gammaTra]^2 * grad2 - Cos[gammaTra]^2;
      (* Return WhenEvent for termination *)
      WhenEvent[vi < 0 // Evaluate, "StopIntegration"]
        /. {x -> x[s], y -> y[s]}
    ]
  ];


(* ::Subsection:: *)
(*Hyperbolic critical terminal point*)


x0[tSol_, gammaTra_] :=
  Module[{p, xMax},
    (* \[PartialD]T/\[PartialD]x *)
    p = Derivative[1, 0][tSol];
    (* Find x == x_0 for which -\[PartialD]T/\[PartialD]x == cot(gamma-tilde) *)
    xMax = 10;
    SeekRoot[p[#, 0] + Cot[gammaTra] &, {0, xMax}]
  ];


(* ::Subsection::Closed:: *)
(*Tracing (old)*)


(* ::Subsubsection:: *)
(*Representative angles*)


(* alpha, gamma and gamma-tilde (tracing gamma) per degree *)
apdRep = 40;
gpdRep = 60;
gpdTraRep = 70;


alphaRep = apdRep * Degree;
gammaRep = gpdRep * Degree;
gammaTraRep = gpdTraRep * Degree;


(* ::Subsubsection:: *)
(*Import known solution*)


tSolRep = Import[
  FString @ "solution/wedge_acute-solution-apd-{apdRep}-gpd-{gpdRep}.txt"
] // Uncompress // First;


meshRep = tSolRep["ElementMesh"];


(* ::Subsubsection:: *)
(*Same contact angle for tracing*)


(* ::Text:: *)
(*Hyperbolic critical terminal point (x_0)*)


x0Same = x0[tSolRep, gammaRep];


(* ::Text:: *)
(*Traced boundaries through (x, y) = (x_0, 0)*)


(* (This is somewhat slow (~2 sec), so compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["same/wedge_acute_same-traced-hyperbolic.txt",
  List @ With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{tSol, gamma, xInit, sMax},
      tSol = tSolRep;
      gamma = gammaRep;
      xInit = x0Same;
      (* Numerical integration *)
      sMax = 3 x0Same;
      NDSolveValue[
        {
          xyTraSystem[tSol, gamma],
          x[0] == xInit,
          y[0] == 0
        }, {x, y}, {s, -sMax, sMax},
        NoExtrapolation
      ] // ReInterpolate
    ]
  ] // Compress
]


(* ::Text:: *)
(*Local T-contour through (x, y) = (x_0, 0)*)


(* (This is somewhat slow (~2 sec), so compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["same/wedge_acute_same-contour-hyperbolic.txt",
  List @ With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{tSol, xInit, sMax},
      tSol = tSolRep;
      xInit = x0Same;
      (* Numerical integration *)
      sMax = 3 x0Same;
      NDSolveValue[
        {
          tContourSystem[tSol],
          x[0] == xInit,
          y[0] == 0
        }, {x, y}, {s, -sMax, sMax},
        NoExtrapolation
      ] // ReInterpolate
    ]
  ] // Compress
]


(* ::Text:: *)
(*Terminal curve*)


(* (This is somewhat slow (~2 sec), so compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["same/wedge_acute_same-terminal.txt",
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[{tSol, xInit, sMax},
      tSol = tSolRep;
      xInit = x0Same;
      (* Numerical integration *)
      sMax = 3 x0Same;
      NDSolveValue[
        {
          viContourSystem[tSol],
          x[0] == xInit,
          y[0] == 0
        }, {x, y}, {s, -sMax, sMax},
        NoExtrapolation
      ] // ReInterpolate
    ]
  ] // Compress
]


(* ::Text:: *)
(*Generic traced boundaries*)


(* (This is slow (~20 sec), so compute once and store.) *)
(* (Delete the file manually to compute from scratch.) *)
ExportIfNotExists["same/wedge_acute_same-traced-general.txt",
  With[{x = \[FormalX], y = \[FormalY], s = \[FormalS]},
    Module[
     {tSol, alpha, gamma, sMax,
      xMax, yMax,
      num, xyInitList, xInit, yInit,
      xyTraWall, xyTraTerm
     },
      tSol = tSolRep;
      alpha = alphaRep;
      gamma = gammaRep;
      (* Numerical integration *)
      sMax = 3 x0Same;
      (* Ranges *)
      xMax = 2 x0Same;
      yMax = xMax Tan[alpha];
      (* Traced boundaries from the lower wall *)
      num = 16;
      xyInitList = Subdivide[{0, 0}, {xMax, -yMax}, num] // Rest;
      xyTraWall =
        Table[
          {xInit, yInit} = xyInit;
          NDSolveValue[
            {
              xyTraSystem[tSol, gamma, -1],
              x[0] == xInit,
              y[0] == yInit,
              xyTraSystemTerm[tSol, gamma]
            }, {x, y}, {s, 0, sMax},
            NoExtrapolation
          ] // ReInterpolate
        , {xyInit, xyInitList}];
      (* Traced boundaries from terminal curve *)
      xyInitList =
        Cases[
          (* Ending points of xyTraWall *)
          Table[
            xy @ DomainEnd[xy] // Through
          , {xy, xyTraWall}],
          (* Reflect those which are on the terminal curve *)
          (* (due to traced boundary hitting terminal curve and terminating) *)
          {x_, y_} /; y < 0
            :>
          {x, -y}
        ];
      xyTraTerm =
        Table[
          {xInit, yInit} = xyInit;
          NDSolveValue[
            {
              xyTraSystem[tSol, gamma, -1],
              x[0] == xInit,
              y[0] == yInit,
              xyTraSystemTerm[tSol, gamma]
            }, {x, y}, {s, 0, sMax},
            NoExtrapolation
          ] // ReInterpolate
        , {xyInit, xyInitList}];
      (* Join lists *)
      Join[xyTraWall, xyTraTerm]
    ]
  ] // Compress
]


(* ::Subsection::Closed:: *)
(*Both branches of a curve (old)*)


(* Adds reflection about x == 0. *)
bothBranches[{x_, y_}] := {{x, y}, {x, -y}};


(* ::Subsection:: *)
(*Geometric regions*)


(* ::Subsubsection:: *)
(*Wedge walls*)


wedgeWalls[alpha_] := {
  HalfLine[{0, 0}, XYPolar[1, alpha]],
  HalfLine[{0, 0}, XYPolar[1, -alpha]]
};


(* ::Subsection:: *)
(*Global styles for plots*)


contStyle = Directive[Thin, Black];
nonStyle = Directive[Opacity[0.7], LightGray];
termStyle = Pink;
wallStyle = Directive[Thick, Purple];


upperStyle = Blue;
lowerStyle = Red;


(* ::Section:: *)
(*Parsing of boundary value problem*)


(* ::Text:: *)
(*Inner workings of SolveLaplaceYoung in LaplaceYoung.wl.*)


Module[
  {
    gamma = 1 Degree,
    mesh = ToElementMesh @ Disk[],
    prWet = True &,
    dummyForTrailingCommas
  },
  With[{x = \[FormalX], y = \[FormalY], t = \[FormalCapitalT]},
    Module[{grad, iGrad, iDiv, k, stateData, propList},
      grad = Grad[#, {x, y}] &;
      iGrad = Inactive[Grad][#, {x, y}] &;
      iDiv = Inactive[Div][#, {x, y}] &;
      k = 1 / Sqrt[1 + # . #] & @ grad @ t[x, y];
      (* Parsing *)
      stateData =
        Quiet[
          NDSolve`ProcessEquations[
            {
              iDiv[-k * IdentityMatrix[2] . iGrad @ t[x, y]] ==
                - t[x, y]
                + NeumannValue[Cos[gamma], prWet[x, y]]
            }
            , t
            , Element[{x, y}, mesh]
          ]
          , {NDSolve`ProcessEquations::femibcnd}
        ];
      propList = {
        "ReactionCoefficients",
        "DiffusionCoefficients",
        Nothing
      };
      Table[
        {
          prop,
          stateData[[1]]["FiniteElementData"]["PDECoefficientData"][prop]
        }
        , {prop, propList}
      ] // TableForm
    ]
  ]
]


(* ::Section:: *)
(*Finite element mesh*)


(* For each value of alpha: *)
Table[
  Module[{mesh, prWet, nElem},
    (* Import mesh *)
    {mesh, prWet} =
      Import[
        FString @ "mesh/wedge_acute-mesh-apd-{apd}.txt"
      ] // Uncompress;
    nElem = Length @ mesh[[2, 1, 1]];
    (* Plot *)
    Show[
      mesh["Wireframe"],
      PlotLabel -> FString @ "{nElem} mesh elements",
      PlotOptions[Axes] // Evaluate
    ]
  ] // Ex @ FString @ "mesh/wedge_acute-mesh-apd-{apd}.png"
, {apd, apdValues}]


(* ::Section:: *)
(*Numerical solution*)


(* ::Subsection:: *)
(*Using built-in nonlinear solver*)


(* ::Subsubsection:: *)
(*Table*)


(* (This is slow (~1 min) from all the calls to Import.) *)
Join @@ Table[
  Table[
    Module[{tSol, messages, timing, tCorner},
      (* Import numerical solution, error messages and abs. time *)
      {tSol, messages, timing} = Import[
        FString @ "solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"
      ] // Uncompress // Part[#, {1, 4, 5}] &;
      (* Filter out uninformative messages and make pretty *)
      messages = DeleteCases[
        messages,
        str_ /; StringContainsQ[str, "StringForm"]
      ];
      messages = StringRiffle[messages, "\n"];
      messages = messages /. {"" -> "\[HappySmiley]"};
      (* Compute corner height *)
      tCorner = tSol[0, 0];
      (* Build row of table *)
      {apd * Degree, gpd * Degree, tCorner, timing, messages}
    ]
  , {gpd, gpdValues}]
, {apd, apdValues}] // TableForm[#,
  TableAlignments -> {Left, Center},
  TableDepth -> 2,
  TableHeadings -> {
    None,
    {"\[Alpha]", "\[Gamma]", "T_corner", "Abs. time", "Messages"}
  },
  TableSpacing -> {3, 3}
] & // Ex["wedge_acute-solution-table.pdf"]


(* ::Subsubsection:: *)
(*Table for theoretical vs computed values*)


(* (This is slow (~6 sec) from all the calls to Import.) *)
Module[
  {
    apd, alpha,
    gamma, k,
    infiniteHeightQ, infiniteSlopeQ,
    theoreticalHeight, theoreticalSlope,
    tSol,
    computedHeight, computedSlope, slopeRelError,
    dummyForTrailingCommas
  },
  (* Representative value of wedge half-angle *)
  apd = apdRep;
  alpha = alphaRep;
  (* For each contact angle *)
  Table[
    (* Theoretical corner height and slope *)
    gamma = gpd Degree;
    k = Sin[alpha] / Cos[gamma];
    infiniteHeightQ = alpha < Pi/2 - gamma;
    infiniteSlopeQ = alpha <= Pi/2 - gamma;
    theoreticalHeight = If[infiniteHeightQ, Infinity, "finite"];
    theoreticalSlope = If[infiniteSlopeQ, Infinity, 1 / Sqrt[k^2 - 1] // N];
    (* Import numerical solution *)
    tSol =
      StringTemplate["solution/wedge_acute-solution-apd-``-gpd-``.txt"][apd, gpd]
        // Import // Uncompress // First;
    (* Computed corner height and slope *)
    computedHeight = tSol[0, 0];
    computedSlope = -Derivative[1, 0][tSol][0, 0];
    slopeRelError = If[infiniteSlopeQ, "N/A", computedSlope / theoreticalSlope - 1];
    (* Build row of table *)
    {
      gamma,
      theoreticalHeight,
      computedHeight,
      theoreticalSlope,
      computedSlope,
      slopeRelError,
      Nothing
    }
    , {gpd, gpdValues}
  ]
    //
      TableForm[#
        , TableHeadings -> {
            None,
            {
              "\[Gamma]",
              "height",
              "(comp.)",
              "slope",
              "(comp.)",
              "(rel. error)",
              Nothing
            }
          }
      ] &
    // Ex["wedge_acute-height-slope-comparison.pdf"]
]


(* ::Subsubsection:: *)
(*Borderline case comparison with asymptotic*)


(* ::Text:: *)
(*From (34) of King et al. (1999) Quart. J. Mech. Appl. Math., 52 (1), 73\[Dash]97:*)
(*T ~ T_0 - 2 sqrt(x/T_0) + y^2/4 sqrt(T_0/x)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tSol, t0,
    tAsy,
    phiValues,
    dummyForTrailingCommas
  },
  (* Borderline case *)
  apd = 40;
  gpd = 90 - apd;
  (* Angles *)
  alpha = apd * Degree;
  gamma = gpd * Degree;
  (* Import numerical solution *)
  tSol =
    Import @ StringTemplate["solution/wedge_acute-solution-apd-``-gpd-``.txt"][apd, gpd]
      // Uncompress // First;
  (* Corner height *)
  t0 = tSol[0, 0];
  (* Asymptotic solution *)
  tAsy[x_, y_] := t0 - 2 Sqrt[x / t0] + y^2/4 Sqrt[t0/x];
  (* Plot discrepancy *)
  phiValues = Subdivide[0, alpha, 4];
  Table[
    Plot[
      {tSol[x, x Tan[phi]], tAsy[x, x Tan[phi]]} // Evaluate
      , {x, 0, 3}
      , ImageSize -> 240
      , PlotLabel -> HoldForm @ "\[Phi]" == phi
      , PlotRange -> Full
      , PlotStyle -> {Automatic, Dashed}
      , PlotOptions[Axes] // Evaluate
    ]
    , {phi, phiValues}
  ]
] // Ex["wedge_acute-borderline-asymptotic-comparison-test.pdf"]


(* ::Section::Closed:: *)
(*Traced boundary plots (old)*)


(* ::Subsection:: *)
(*Same contact angle for tracing*)


(* ::Subsubsection:: *)
(*Hyperbolic traced boundaries*)


Module[
 {alpha, gamma,
  xMax, yMax,
  xyTerm, xyCont, xyTra
 },
  (* Constants *)
  alpha = alphaRep;
  gamma = gammaRep;
  (* Plot range *)
  xMax = 2 x0Same;
  yMax = xMax Tan[alpha];
  (* Import curves *)
  xyTerm = Import @ "same/wedge_acute_same-terminal.txt" // Uncompress;
  xyCont = Import @ "same/wedge_acute_same-contour-hyperbolic.txt" // Uncompress;
  xyTra = Import @ "same/wedge_acute_same-traced-hyperbolic.txt" // Uncompress;
  (* Plot *)
  Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}],
    (* Walls *)
    Graphics @ {wallStyle, wedgeWalls[alpha]},
    (* Non-viable domain and terminal curve *)
    Table[
      ParametricPlot[
        {
          Way[xy[[1]][s], 2 xMax, p],
          xy[[2]][s]
        },
        {s, DomainStart[xy], DomainEnd[xy]},
        {p, 0, 1},
        PlotStyle -> nonStyle,
        BoundaryStyle -> termStyle
      ]
    , {xy, {xyTerm}}],
    (* Local T-contour *)
    Table[
      ParametricPlot[
        xy[s] // Through,
        {s, DomainStart[xy], DomainEnd[xy]},
        PlotStyle -> contStyle
      ]
    , {xy, xyCont}],
    (* Traced boundaries *)
    Table[
      ParametricPlot[
        xy[s] // Through // bothBranches // Evaluate,
        {s, DomainStart[xy], DomainEnd[xy]},
        PlotStyle -> {upperStyle, lowerStyle}
      ]
    , {xy, xyTra}]
  ]
] // Ex @ FString @ StringJoin[
  "same/wedge_acute_same-traced-hyperbolic",
  "-apd-{apdRep}-gpd-{gpdRep}-gpdt-{gpdRep}.pdf"
]


(* ::Subsubsection:: *)
(*Generic traced boundaries*)


Module[
 {alpha, gamma,
  xMax, yMax,
  xyTerm, xyCont, xyTra
 },
  (* Constants *)
  alpha = alphaRep;
  gamma = gammaRep;
  (* Plot range *)
  xMax = 2 x0Same;
  yMax = xMax Tan[alpha];
  (* Import curves *)
  xyTerm = Import @ "same/wedge_acute_same-terminal.txt" // Uncompress;
  xyCont = Import @ "same/wedge_acute_same-contour-hyperbolic.txt" // Uncompress;
  xyTra = Import @ "same/wedge_acute_same-traced-general.txt" // Uncompress;
  (* Plot *)
  Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}],
    (* Walls *)
    Graphics @ {wallStyle, wedgeWalls[alpha]},
    (* Non-viable domain and terminal curve *)
    Table[
      ParametricPlot[
        {
          Way[xy[[1]][s], 2 xMax, p],
          xy[[2]][s]
        },
        {s, DomainStart[xy], DomainEnd[xy]},
        {p, 0, 1},
        PlotStyle -> nonStyle,
        BoundaryStyle -> termStyle
      ]
    , {xy, {xyTerm}}],
    (* Traced boundaries *)
    Table[
      ParametricPlot[
        xy[s] // Through // bothBranches // Evaluate,
        {s, DomainStart[xy], DomainEnd[xy]},
        PlotStyle -> {upperStyle, lowerStyle}
      ]
    , {xy, xyTra}]
  ]
] // Ex @ FString @ StringJoin[
  "same/wedge_acute_same-traced-general",
  "-apd-{apdRep}-gpd-{gpdRep}-gpdt-{gpdRep}.pdf"
]


(* ::Section:: *)
(*Figure: finite element mesh (wedge_acute-mesh-*)*)


(* ::Subsection:: *)
(*Full*)


(* (This is slow (~6 sec).) *)
Module[
  {
    apd, alpha,
    mesh, meshWireframe,
    wallDistanceSimplified, polygonSimplified,
    posSimplified, meshWireframeSimplified,
    xMin, xMax, yMin, yMax, rMax,
    tickTextStyle, axesTextStyle,
    phiPattern,
    plotFull,
    dummyForTrailingCommas
  },
  (* Contact angle *)
  apd = 40;
  alpha = apd * Degree;
  (* Import mesh *)
  mesh = Import @ FString["mesh/wedge_acute-mesh-apd-{apd}.txt"] // Uncompress // First;
  meshWireframe = mesh["Wireframe"];
  {{xMin, xMax}, {yMin, yMax}} = mesh["Bounds"];
  rMax = xMax;
  (*
    Simplify mesh to reduce file size
    by replacing detailed portions with black polygons
  *)
  wallDistanceSimplified = 0.1;
  polygonSimplified = Graphics @ Polygon @ {
    {wallDistanceSimplified / Sin[alpha], 0},
    XYPolar[rMax, alpha] + XYPolar[wallDistanceSimplified, alpha - Pi/2],
    XYPolar[rMax, alpha],
    {0, 0},
    XYPolar[rMax, -alpha],
    XYPolar[rMax, -alpha] + XYPolar[wallDistanceSimplified, -alpha + Pi/2],
    Nothing
  };
  posSimplified =
    MeshWireframePositions[mesh,
      {x_, y_} /;
        x Sin[alpha] - y Cos[alpha] > wallDistanceSimplified
          &&
        x Sin[alpha] + y Cos[alpha] > wallDistanceSimplified
    ];
  meshWireframeSimplified = mesh["Wireframe"[posSimplified]];
  (* Plot full mesh *)
  tickTextStyle = Style[#, 13] & @* LaTeXStyle;
  axesTextStyle = Style[#, 16] & @* LaTeXStyle;
  phiPattern = _Integer Degree | 0;
  plotFull = Show[
    PolarPlot[rMax, {phi, -alpha, alpha}
      , PlotStyle -> None
      , PolarAxes -> {True, True}
      , PolarAxesOrigin -> {-alpha, rMax}
      , PolarTicks -> {Range[-alpha, alpha, 10 Degree], Automatic}
    ]
      (* Angular PlotRange *)
      /. {Circle[seq__, {0, 2 Pi}] :> Circle[seq, {-alpha, alpha}]}
      (* Tweak radial labels *)
      /. {Text[Style[r : Except[phiPattern], {}], coords_, offset_] :>
        Text[r, coords, offset + If[r == 10, {1.5, 0}, {0.75, 0}]]
      }
      (* Tweak azimuthal labels *)
      /. {Text[Style[phi : phiPattern, {}], coords_, offset_, opts___] :>
        Text[
          SeparatedRow["VeryThin"][phi / Degree, Magnify["\[Degree]", 1.2]],
          coords,
          offset + Which[
            phi == -alpha, {-0.2, -0.7},
            phi < 0, Abs[phi]/alpha {0, -0.7},
            phi == 0, {-0.2, 0},
            True, {0, 0}
          ]
          , opts
        ]
      }
      (* Font for labels *)
      /. {Text[str_, seq__] :> Text[str // tickTextStyle, seq]}
    ,
    meshWireframeSimplified, polygonSimplified,
    (* Manual coordinate labels *)
    Graphics @ {
      Text[
        Italicise["r"] // axesTextStyle
        , XYPolar[rMax / 2, -alpha]
        , {5.5, 2.2}
      ],
      Text[
        "\[Phi]" // axesTextStyle // Margined[{{0, 1}, {0, 0}}]
        , {rMax, 0}
        , {-8, 0.1}
      ],
      {}
    },
    {}
    , ImageSize -> 360
    , PlotRange -> All
    , PlotRangeClipping -> False
  ];
  plotFull // Ex["wedge_acute-mesh-full.pdf"]
]


(* ::Subsection:: *)
(*Detail*)


Module[
  {
    apd, alpha,
    mesh, meshWireframe,
    rMax, xMax, yMax,
    tickTextStyle, axesTextStyle,
    plotDetail,
    dummyForTrailingCommas
  },
  (* Contact angle *)
  apd = 40;
  alpha = apd * Degree;
  (* Import mesh *)
  mesh = Import @ FString["mesh/wedge_acute-mesh-apd-{apd}.txt"] // Uncompress // First;
  meshWireframe = mesh["Wireframe"];
  rMax = 0.25;
  {xMax, yMax} = XYPolar[rMax, alpha];
  (* Plot mesh detail *)
  tickTextStyle = Style[#, 15] & @* LaTeXStyle;
  axesTextStyle = Style[#, 18] & @* LaTeXStyle;
  plotDetail = Show[
    PolarPlot[rMax, {phi, -alpha, alpha}
      , PlotStyle -> None
      , PolarAxes -> {False, True}
      , PolarAxesOrigin -> {-alpha, rMax}
      , PolarTicks -> List @
          Table[
            {r, N[r], {0, 0.04 rMax}}
            , {r, FindDivisions[{0, rMax}, 4]}
          ]
    ]
      (* Tweak radial labels *)
      /. {Text[Style[r_, {}], coords_, offset_, opts__] :>
        Text[r, coords, offset + {If[r == 0, 0.8, 1.4], 0.1}, opts]
      }
      (* Font for labels *)
      /. {Text[str_, seq__] :> Text[str // tickTextStyle, seq]}
    ,
    meshWireframe,
    (* Manual coordinate labels *)
    Graphics @ {
      Text[
        Italicise["r"] // axesTextStyle
        , XYPolar[rMax / 2, -alpha]
        , {3, 2.5}
      ],
      {}
    },
    {}
    , ImageSize -> 180
    , PlotRange -> {{0, xMax}, {-yMax, yMax}}
    , PlotRangePadding -> {0.2, {0.3, 0.2}} rMax
  ];
  plotDetail
] // Ex["wedge_acute-mesh-detail.pdf"]


(* ::Section:: *)
(*Figure: borderline case comparison with asymptotic  (wedge_acute-borderline-asymptotic-comparison-*)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical, t0,
    tAsymptotic,
    styleCommon, styleNumerical, styleAsymptotic,
    phiValues, phiCases, case,
    dummyForTrailingCommas
  },
  (* Borderline case *)
  apd = 40;
  gpd = 90 - apd;
  (* Angles *)
  alpha = apd * Degree;
  gamma = gpd * Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Corner height *)
  t0 = tNumerical[0, 0];
  (* Asymptotic solution *)
  tAsymptotic[x_, y_] := t0 - 2 Sqrt[x / t0] + y^2/4 Sqrt[t0/x];
  (* Plot styles *)
  styleCommon = Directive[Black, GeneralStyle["DefaultThick"]];
  styleNumerical = styleCommon;
  styleAsymptotic = Directive[styleCommon, GeneralStyle["Dashed"]];
  (* Plot discrepancy *)
  phiValues = {0, alpha};
  phiCases = {"symmetry", "wall"};
  Table[
    case = AssociationThread[phiValues -> phiCases][phi];
    Plot[
      {tNumerical[x, x Tan[phi]], tAsymptotic[x, x Tan[phi]]} // Evaluate
      , {x, 0, 3}
      , AxesLabel -> Italicise /@ {"x", "T"}
      , ImageSize -> 270
      , LabelStyle -> LatinModernLabelStyle[14]
      , PlotRange -> {0, All}
      , PlotStyle -> {styleNumerical, styleAsymptotic}
      , PlotOptions[Axes] // Evaluate
    ] // Ex @ FString["wedge_acute-borderline-asymptotic-comparison-{case}.pdf"]
    , {phi, phiValues}
  ]
]


(* ::Section:: *)
(*Figure: known solution  (wedge_acute-solution-*)*)


Module[
  {
    apd, alpha,
    gpdValues, gpdCases,
    rMax, xArcCorner, yArcCorner,
    case,
    tNumerical,
    x, y, mesh,
    dummyForTrailingCommas
  },
  (* Wedge half-angle *)
  apd = 40;
  alpha = apd * Degree;
  (* Contact angles *)
  gpdValues = {50, 55};
  gpdCases = {"borderline", "moderate"};
  (* Plot range *)
  rMax = 4;
  {xArcCorner, yArcCorner} = XYPolar[rMax, alpha];
  (* Make plots *)
  Table[
    case = AssociationThread[gpdValues -> gpdCases][gpd];
    (* Import numerical solution *)
    tNumerical =
      Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
        // Uncompress // First;
    (* Plot *)
    mesh = tNumerical["ElementMesh"];
    Show[
      (* Numerical solution *)
      Plot3D[
        tNumerical[x, y], Element[{x, y}, mesh]
        , AxesEdge -> {{-1, -1}, {+1, -1}, {-1, -1}}
        , AxesLabel -> {
            Italicise["x"],
            Italicise["y"],
            Italicise["T"] // Margined[{{0, 0.2}, {0, 1}} ImageSizeCentimetre]
          }
        , Boxed -> {Back, Bottom, Left}
        , BoxRatios -> Automatic
        , Filling -> 0
        , FillingStyle -> BoundaryTracingStyle["Solution3D"]
        , LabelStyle -> LatinModernLabelStyle[11]
        , Lighting -> GeneralStyle["AmbientLighting"]
        , PlotPoints -> 20
        , PlotRange -> {0, 2.2}
        , PlotStyle -> BoundaryTracingStyle["Solution3D"]
        , TicksStyle -> ConstantArray[LatinModernLabelStyle[8], 3]
        , RegionFunction -> Function[{x, y},
            RPolar[x, y] < rMax
          ]
        , ViewPoint -> {1.5, -2.5, 1.4}
      ],
      (* Manually drawn solution base edges *)
      Graphics3D @ {
        Line @ {{0, 0, 0}, {xArcCorner, -yArcCorner, 0}},
        {}
      },
      ParametricPlot3D[
        rMax {Cos[phi], Sin[phi], 0}
        , {phi, -alpha, alpha}
        , PlotStyle -> Directive[Black, Thickness[0]]
      ],
      (* Manually drawn solution vertical edges *)
      Graphics3D @ {
        Line @ {
          {0, 0, tNumerical[0, 0]},
          {0, 0, 0}
        },
        Line @ {
          {xArcCorner, -yArcCorner, tNumerical[xArcCorner, -yArcCorner]},
          {xArcCorner, -yArcCorner, 0}
        },
        {}
      },
      {}
      , ImageSize -> 6.5 ImageSizeCentimetre
    ]
      // Ex[FString["wedge_acute-solution-{case}.png"]
        , Background -> None
        , ImageResolution -> 4 BasicImageResolution
      ]
    , {gpd, gpdValues}
  ]
]


(* ::Section:: *)
(*Figure: generic traced boundaries  (wedge_acute-traced-boundaries)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical, xCritical,
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    derList, p, q, grad2, f, vi,
    rStartList,
    xyStartList, xyTracedList,
    plot, legend,
    legendLabelStyle,
    legendCurves, legendRegions,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {40, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMax = Ceiling[1.5 xCritical, 0.2];
  yMax = xMax Tan[alpha];
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* Traced boundaries (upper branch) *)
  (*
    NOTE: wall integration doesn't work
    because the traced boundary ventures slightly outside the wedge domain
    and Mathematica complains about Indeterminate values.
    So the wall boundaries are drawn manually at the end.
  *)
  rStartList = Subdivide[rMax, 12] // Rest;
  xyStartList = Table[XYPolar[r, -alpha], {r, rStartList}];
  xyTracedList =
    Table[
      ContactTracedBoundary[derList][xyStart, 0, {0, 2 rMax}
        , -1, 1
        , 0, Function[{x, y}, x > xMaxMore]
      ]
      , {xyStart, xyStartList}
    ];
  (* Plot *)
  plot = Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}
      , ImageSize -> 270
      , LabelStyle -> LatinModernLabelStyle[14]
    ],
    (* Wedge walls *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      HalfLine[{0, 0}, XYPolar[1, alpha]],
      HalfLine[{0, 0}, XYPolar[1, -alpha]],
      {}
    },
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0
      , {x, xCritical, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 7
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Traced boundaries *)
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // IncludeYReflection
          // Evaluate
        , {s, DomainStart[xy], DomainEnd[xy]}
        , PlotStyle -> BoundaryTracingStyle["Traced"]
      ]
      , {xy, xyTracedList}
    ],
    (* Traced boundary walls *)
    ParametricPlot[
      {x, x Tan[alpha]}
        // IncludeYReflection
        // Evaluate
      , {x, 0, xMaxMore}
      , PlotPoints -> 2
      , PlotStyle -> BoundaryTracingStyle["Traced"]
    ],
    {}
  ];
  (* Legend *)
  legendLabelStyle = LatinModernLabelStyle[14];
  legendCurves =
    CurveLegend[
      BoundaryTracingStyle /@ {"Traced", "Terminal"},
      {"traced boundary", "terminal curve"}
      , LabelStyle -> legendLabelStyle
    ];
  legendRegions =
    RegionLegend[
      BoundaryTracingStyle /@ {"NonViable"},
      {"non\[Hyphen]viable domain"}
      , LabelStyle -> legendLabelStyle
    ];
  (* Combined *)
  Grid[
    {{
      plot,
      Column[Join[legendCurves, legendRegions], Spacings -> -0.5]
    }}
    , Spacings -> 2
  ]
] // Ex["wedge_acute-traced-boundaries.pdf"]


(* ::Section:: *)
(*Figure: traced boundaries, patched (wedge_acute-traced-boundaries-patched)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical, xCritical,
    derList, p, q, grad2, f, vi,
    xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    normalVectorLength,
    commonPlot,
    cornerList, numPlots,
    sRange,
    xyTracedUpperList, xyTracedLowerList,
    numCorners,
    sIntersectionLower, sIntersectionUpper,
    xTracedForNormal, yTracedForNormal, sNormal, xyNormal, normalVector,
    xyLower, xyUpper,
    textStyle,
    plotList,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {40, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = derList = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMax = Ceiling[1.5 xCritical, 0.2];
  yMax = xMax Tan[alpha];
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  normalVectorLength = 1/3 xMax;
  (* Location of corners (increasing in y-coordinate) *)
  cornerList = Association[
    1 -> {{0.45, 0.2}},
    2 -> {{0.4, 0}, {0.45, 0.2}},
    3 -> {{0.6, -0.5}, {0.4, 0.25}, {0.5, 0.4}},
    4 -> {{0.5, -0.38}, {0.35, 0}, {0.4, 0.25}, {0.5, 0.4}},
    Nothing
  ];
  numPlots = Length[cornerList];
  (* Compute traced boundaries therethrough *)
  sRange = 2 rMax {-1, 1};
  Table[
    (* Upper branch *)
    xyTracedUpperList[n] =
      Table[
        ContactTracedBoundary[derList][corner, 0, sRange
          , -1, 1
          , 0, Function[{x, y}, x > xMaxMore]
        ]
        , {corner, cornerList[n]}
      ];
    xyTracedLowerList[n] =
      Table[
        ContactTracedBoundary[derList][corner, 0, sRange
          , 1, -1
          , 0, Function[{x, y}, x > xMaxMore]
        ]
        , {corner, cornerList[n]}
      ];
    (* Intersections *)
    numCorners = Length @ cornerList[n];
    Table[
      (* Starting arc length for lower branch *)
      sIntersectionLower[n][k] =
        If[k == 1
          ,
          xyLower = xyTracedLowerList[n][[k]];
          DomainEnd[xyLower]
          ,
          xyLower = xyTracedLowerList[n][[k]];
          xyUpper = xyTracedUpperList[n][[k - 1]];
          SeekParametricIntersection[xyLower, xyUpper] // First
        ];
      (* Ending arc length for upper branch *)
      sIntersectionUpper[n][k] =
        If[k == numCorners
          ,
          xyUpper = xyTracedUpperList[n][[k]];
          DomainEnd[xyUpper]
          ,
          xyUpper = xyTracedUpperList[n][[k]];
          xyLower = xyTracedLowerList[n][[k + 1]];
          SeekParametricIntersection[xyUpper, xyLower] // First
        ];
      , {k, numCorners}
    ];
    (* Normal vector *)
    {xTracedForNormal, yTracedForNormal} = xyTracedUpperList[n][[1]];
    sNormal = Way[0, sIntersectionUpper[n][1], If[n == 1, 1/4, 1/2]];
    xyNormal[n] = {xTracedForNormal, yTracedForNormal}[sNormal] // Through;
    normalVector[n] = (
      {xTracedForNormal', yTracedForNormal'}[sNormal]
        // Through
        // Normalize
        // Cross
    );
    , {n, numPlots}
  ];
  (* Common plot *)
  commonPlot = Show[
    EmptyFrame[{0, xMax}, {-yMax, yMax}
      , Frame -> None
      , ImageSize -> 180
    ],
    {}
  ];
  (* Plots *)
  textStyle = Style[#, 24] & @* LaTeXStyle;
  plotList = Table[
    Show[
      (* Common plot *)
      commonPlot,
      (* Traced boundaries (upper) *)
      Table[
        ParametricPlot[
          xy[s]
            // Through
            // Evaluate
          , {s, DomainStart[xy], DomainEnd[xy]}
          , PlotStyle -> BoundaryTracingStyle["Background"]
        ]
        , {xy, xyTracedUpperList[n]}
      ],
      (* Traced boundaries (background) *)
      Table[
        ParametricPlot[
          xy[s]
            // Through
            // Evaluate
          , {s, DomainStart[xy], DomainEnd[xy]}
          , PlotStyle -> BoundaryTracingStyle["Background"]
        ]
        , {xy, xyTracedLowerList[n]}
      ],
      (* Normal vector *)
      Graphics @ {Directive[GeneralStyle["Thick"], Arrowheads[Large]],
        Arrow @ {
          xyNormal[n],
          xyNormal[n] + normalVectorLength * normalVector[n]
        }
      },
      Graphics @ {
        Text[
          Embolden["n"] // textStyle
          , xyNormal[n] + normalVectorLength * normalVector[n]
          , {1., -0.7}
        ]
      },
      (* Traced boundaries (patched) *)
      numCorners = Length @ cornerList[n];
      Table[
        {
          (* Lower branch portions *)
          ParametricPlot[
            xyTracedLowerList[n][[k]][s]
              // Through
              // Evaluate
            , {s, sIntersectionLower[n][k], 0}
            , PlotStyle -> BoundaryTracingStyle["Traced"]
          ],
          (* Upper branch portions *)
          ParametricPlot[
            xyTracedUpperList[n][[k]][s]
              // Through
              // Evaluate
            , {s, 0, sIntersectionUpper[n][k]}
            , PlotStyle -> BoundaryTracingStyle["Traced"]
          ],
          {}
        }
        , {k, numCorners}
      ],
      {}
    ]
  , {n, numPlots}];
  Grid @ {plotList}
] // Ex["wedge_acute-traced-boundaries-patched.pdf"]


(* ::Section:: *)
(*Figure: terminal-points (wedge_acute-terminal-points)*)


Module[
  {
    apd, gpd,
    alpha, gamma,
    tNumerical, xCritical,
    p, q, grad2, f, vi,
    xMin, xMax, yMax, rMax,
    more, xMaxMore, yMaxMore,
    xyStartList, xyContourList,
    xyTerminal,
    xyContourOrdinary, sOrdinary, xyOrdinary,
    textStyle, textStyleBracket, textVerticalShift,
    plot,
    legendLabelStyle,
    legendCurves, legendRegions,
    dummyForTrailingCommas
  },
  (* Angular parameters *)
  {apd, gpd} = {40, 60};
  {alpha, gamma} = {apd, gpd} Degree;
  (* Import numerical solution *)
  tNumerical =
    Import @ FString["solution/wedge_acute-solution-apd-{apd}-gpd-{gpd}.txt"]
      // Uncompress // First;
  (* Derivative list for boundary tracing *)
  {p, q, grad2, f, vi} = ContactDerivativeList[tNumerical, gamma];
  (* Critical terminal point x_0 *)
  xCritical = x0[tNumerical, gamma];
  (* Plot range *)
  xMin = xCritical / 2;
  xMax = 2 xCritical;
  yMax = xMax Tan[alpha];
  rMax = RPolar[xMax, yMax];
  (* Plot range but more *)
  more = 0.05;
  xMaxMore = xMax + more;
  yMaxMore = yMax + more;
  (* T-contours *)
  xyStartList = Table[{x, 0}, {x, xCritical + 0.07 Range[0, 3]}];
  xyContourList =
    Table[
      ContourByArcLength[tNumerical][xyStart, 0, 2 rMax {-1, 1}, 1, #1 > xMaxMore &]
      , {xyStart, xyStartList}
    ];
  (* Terminal curve *)
  xyTerminal =
    ContourByArcLength[vi][{xCritical, 0}, 0, 2 rMax {-1, 1}, 1, #1 > xMaxMore &];
  (* Orindary terminal point to be shown *)
  xyContourOrdinary = xyContourList[[2]];
  sOrdinary =
    First @ SeekParametricIntersection[
      xyContourOrdinary, xyTerminal,
      {0, DomainEnd[xyContourOrdinary]}, {0, DomainEnd[xyTerminal]}
    ];
  xyOrdinary = xyContourOrdinary[sOrdinary] // Through;
  (* Text style *)
  textStyle = Style[#, 18] & @* LaTeXStyle;
  textStyleBracket = Style[#, Larger] &;
  textVerticalShift = -0.25;
  (* Plot *)
  plot = Show[
    EmptyFrame[{xMin, xMax}, {-yMax, yMax}
      , Frame -> None
      , ImageSize -> 210
    ],
    (* Contours *)
    Table[
      ParametricPlot[
        xy[s]
          // Through
          // Evaluate
        , {s, DomainStart[xy], DomainEnd[xy]}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Background"]
      ]
      , {xy, xyContourList}
    ],
    (* Non-viable domain *)
    RegionPlot[
      vi[x, y] < 0
      , {x, 0, xMaxMore}
      , {y, -yMaxMore, yMaxMore}
      , BoundaryStyle -> BoundaryTracingStyle["Terminal"]
      , PlotPoints -> 7
      , PlotStyle -> BoundaryTracingStyle["NonViable"]
    ],
    (* Critical terminal point (x_0, 0) *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ {xCritical, 0}
    },
    Graphics @ {
      Text[
        Row @ {
          "(" // textStyleBracket,
          "\[NegativeVeryThinSpace]",
          Subscript[Italicise["x"], 0],
          ",\[ThinSpace]",
          0,
          ")" // textStyleBracket
        },
        {xCritical, 0},
        {1.5, textVerticalShift}
      ] // textStyle,
      Text[
        "critical",
        {xCritical, 0},
        {-1.4, textVerticalShift}
      ] // textStyle,
      {}
    },
    (* Ordinary terminal point *)
    Graphics @ {
      GeneralStyle["Point"],
      Point @ xyOrdinary
    },
    Graphics @ {
      Text[
        "ordinary",
        xyOrdinary,
        {-1.4, textVerticalShift}
      ] // textStyle,
      {}
    },
    {}
  ];
  (* Legend *)
  legendLabelStyle = LatinModernLabelStyle[16];
  legendCurves =
    CurveLegend[
      BoundaryTracingStyle /@ {"Background", "Terminal"},
      {Row @ {Italicise["T"], "\[Hyphen]contour"}, "terminal curve"}
      , LabelStyle -> legendLabelStyle
    ];
  legendRegions =
    RegionLegend[
      BoundaryTracingStyle /@ {"NonViable"},
      {"non\[Hyphen]viable domain"}
      , LabelStyle -> legendLabelStyle
    ];
  (* Combined *)
  Grid[
    {{
      plot,
      Column[Join[legendCurves, legendRegions], Spacings -> -0.5]
    }}
    , Spacings -> 3
  ]
] // Ex["wedge_acute-terminal-points.pdf"]
