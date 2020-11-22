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
    textStyle,
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
  textStyle = Style[#, 15] & @* LaTeXStyle;
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
          phi,
          coords,
          offset + Which[
            phi == -alpha, {-0.2, -0.7},
            phi < 0, Abs[phi]/alpha {0, -0.7},
            True, {0, 0}
          ]
          , opts
        ]
      }
      (* Font for labels *)
      /. {Text[str_, seq__] :> Text[str // textStyle, seq]}
    ,
    meshWireframeSimplified, polygonSimplified,
    (* Manual coordinate labels *)
    Graphics @ {
      Text[
        Italicise["r"] // textStyle
        , XYPolar[rMax / 2, -alpha]
        , 2.7 {2, 1}
      ],
      Text[
        "\[Phi]" // textStyle // Margined[{{0, 1}, {0, 0}}]
        , {rMax, 0}
        , {-8, 0}
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
