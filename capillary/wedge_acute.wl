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


(* ::Subsection:: *)
(*System of ODEs for T-contour*)


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


(* ::Subsection:: *)
(*System of ODEs for terminal curve \[CapitalPhi] = 0*)


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


(* ::Subsection:: *)
(*System of ODEs for tracing*)


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


(* ::Subsection:: *)
(*Tracing*)


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
      ]
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
      ]
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
          ]
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
          ]
        , {xyInit, xyInitList}];
      (* Join lists *)
      Join[xyTraWall, xyTraTerm]
    ]
  ] // Compress
]


(* ::Subsection:: *)
(*Both branches of a curve*)


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


(* ::Section:: *)
(*Traced boundary plots*)


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
  "same/wedge_acute-same-traced-hyperbolic",
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
