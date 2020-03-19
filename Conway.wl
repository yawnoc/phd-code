(* ::Package:: *)

(* ::Text:: *)
(*Conway.wl (for Wolfram Mathematica)*)


(* ::Text:: *)
(*A Wolfram Language Package (.wl) with handy functions.*)
(*(C) 2019 Conway*)
(*ABSOLUTELY NO WARRANTY, i.e. "GOD SAVE YOU"*)


(* ::Section:: *)
(*Improve workflow*)


(* ::Text:: *)
(*Disable "smart" quotes and annoying floating elements*)
(*(which pop up when I don't need them),*)
(*and make alternative input auto replacements*)
(*(to avoid unreadable \[...] representations in text editors).*)


SetOptions[$FrontEndSession,
  AutoQuoteCharacters -> {},
  CodeAssistOptions -> {"FloatingElementEnable" -> False},
  InputAutoReplacements -> {
    "eee" -> "==",
    "ggg" -> ">=",
    "lll" -> "<=",
    "nnn" -> "!=",
    "ddd" -> ":>",
    "rrr" -> "->"
  }
];


(* ::Section:: *)
(*Package definitions*)


(* ::Subsection:: *)
(*Start of package*)


BeginPackage["Conway`"];


(* ::Subsection:: *)
(*Clear existing definitions if any*)


(* ::Text:: *)
(*https://mathematica.stackexchange.com/a/128790*)
(*https://mathematica.stackexchange.com/a/57506*)


Unprotect["Conway`*"];
ClearAll["Conway`*"];
ClearAll["Conway`*`*"];


(* ::Subsection:: *)
(*Mention non-private symbols*)


{
  BoxedLabel,
  DefaultColours,
  DeleteNearbyPoints,
  DomainEnd,
  DomainStart,
  EmptyAxes,
  EmptyFrame,
  Ex,
  ExportIfNotExists,
  FString,
  Italicised,
  ModuleSymbol,
  NoExtrapolation,
  PlotOptions,
  PreciseOptions,
  PrettyString,
  ReInterpolate,
  SeekFirstRootBisection,
  SeekRoot,
  SeekRootBisection,
  SortByPhi,
  ToName,
  UniformRange,
  Way
};


(* ::Subsection:: *)
(*Private scope*)


(* ::Subsubsection:: *)
(*Start of private scope*)


Begin["`Private`"];


(* ::Subsubsection:: *)
(*BoxedLabel*)


BoxedLabel::usage = (
  "BoxedLabel[expr, opts]\n"
  <> "Returns framed version of the expression expr.\n"
  <> "To be used for PlotLabel."
);


BoxedLabel[expr_, opts : OptionsPattern[Framed]] :=
  Framed[expr, ImageMargins -> {{0, 0}, {10, 0}}, opts];


(* ::Subsubsection:: *)
(*DefaultColours*)


DefaultColours::usage = (
  "DefaultColours\n"
  <> "Returns the default colours for Plot since Version 10."
);


DefaultColours = ColorData[97, "ColorList"];


(* ::Subsubsection:: *)
(*DeleteNearbyPoints*)


DeleteNearbyPoints::usage = (
  "DeleteNearbyPoints[sep][list]\n"
  <> "Removes duplicates from the list of points list, "
  <> "where points with separation less than sep are considered duplicate."
);


DeleteNearbyPoints[sep_?NumericQ][list_List] :=
  DeleteDuplicates[list, Norm[#2 - #1] < sep &];


(* ::Subsubsection:: *)
(*DomainEnd*)


DomainEnd::usage = (
  "DomainEnd[iFun] or DomainEnd[{iFun, ...}]\n"
  <> "Returns the end of the domain for the first coordinate "
  <> "of the interpolating function iFun."
);


DomainEnd[iFun_InterpolatingFunction] := iFun["Domain"][[1, 2]];


DomainEnd[{iFun_, ___}] := DomainEnd[iFun];


(* ::Subsubsection:: *)
(*DomainStart*)


DomainStart::usage = (
  "DomainStart[iFun] or DomainStart[{iFun, ...}]\n"
  <> "Returns the start of the domain for the first coordinate "
  <> "of the interpolating function iFun."
);


DomainStart[iFun_InterpolatingFunction] := iFun["Domain"][[1, 1]];


DomainStart[{iFun_, ___}] := DomainStart[iFun];


(* ::Subsubsection:: *)
(*EmptyAxes*)


EmptyAxes::usage = (
  "EmptyAxes[{xMin, xMax}, {yMin, yMax}, opts]\n"
  <> "Returns empty Plot with default options PlotOptions[Axes], "
  <> "which may be overridden by options opts."
);


EmptyAxes[
  {xMin_?NumericQ, xMax_?NumericQ},
  {yMin_?NumericQ, yMax_?NumericQ},
  opts : OptionsPattern[Plot]
] :=
  Plot[Indeterminate, {x, xMin, xMax},
    opts,
    PlotRange -> {{xMin, xMax}, {yMin, yMax}},
    PlotOptions[Axes] // Evaluate
  ];


(* ::Subsubsection:: *)
(*EmptyFrame*)


EmptyFrame::usage = (
  "EmptyFrame[{xMin, xMax}, {yMin, yMax}, opts]\n"
  <> "Returns empty RegionPlot with default options PlotOptions[Frame], "
  <> "which may be overridden by options opts."
);


EmptyFrame[
  {xMin_?NumericQ, xMax_?NumericQ},
  {yMin_?NumericQ, yMax_?NumericQ},
  opts : OptionsPattern[RegionPlot]
] :=
  RegionPlot[False, {x, xMin, xMax}, {y, yMin, yMax},
    opts,
    PlotOptions[Frame] // Evaluate
  ];


(* ::Subsubsection:: *)
(*Ex*)


Ex::usage = (
  "Ex[dest, opts]\n"
  <> "Operator form of Export, equivalent to "
  <> "Export[dest, #, opts] &"
);


Ex[dest_, opts : OptionsPattern[Export]] := Export[dest, #, opts] &;


(* ::Subsubsection:: *)
(*ExportIfNotExists*)


ExportIfNotExists::usage = (
  "ExportIfNotExists[dest, expr, opts]\n"
  <> "Exports the expression expr to the destination dest "
  <> "if the destination file does not already exist."
);


ExportIfNotExists[dest_, expr_, opts : OptionsPattern[Export]] :=
  If[Not @ FileExistsQ[dest],
    expr // Ex[dest, opts]
  ];


SetAttributes[ExportIfNotExists, HoldAll];


(* ::Subsubsection:: *)
(*FString*)


FString::usage = (
  "FString[s]\n"
  <> "Formats string s as per Python 3.6's f-strings: "
  <> "contents of curly brackets are evaluated using ToExpression "
  <> "(use doubled curly brackets for literal brackets).\n"
  <> "Referencing With variables will NOT work.\n"
  <> "Avoid referencing free variables: this will be SLOW.\n"
  <> "Automatically threads over lists.\n"
  <> "See PEP 498 -- Literal String Interpolation: "
  <> "https://www.python.org/dev/peps/pep-0498/"
);


FString[s_String] :=
  s // StringReplace @ {
    (* Doubled curly brackets for literal curly brackets *)
    "{{" -> "{",
    "}}" -> "}",
    (* Evaluate contents of single curly brackets *)
    StringExpression[
      "{",
      contents : Shortest[___],
      "}"
    ] :> ToString[
      ToExpression[contents]
        (* Assume free symbols are Module variables *)
        // ReplaceAll[
          x_Symbol :> ModuleSymbol[SymbolName[x], $ModuleNumber]
        ]
        (* Decrement module number until not free *)
        // ReplaceRepeated[#,
          x_Symbol :> ModuleSymbol[
            (* Name part of Module symbol is everything before $ *)
            SymbolName[x] // StringReplace["$" ~~ __ -> ""],
            (* Number part of Module symbol is everything after $ *)
            SymbolName[x]
              // StringReplace[__ ~~ "$" -> ""]
              // ToExpression
              // Max[1, # - 1] &
              (* (stop when $ModuleNumber == 1) *)
          ],
          MaxIterations -> Infinity
        ] &
        // ReplaceAll[
          x_Symbol :> Symbol[
            SymbolName[x] // StringReplace["$" ~~ __ -> ""]
          ]
        ]
    ]
  };


SetAttributes[FString, Listable];


(* ::Subsubsection:: *)
(*Italicised*)


Italicised::usage = (
  "Italicised[str]\n"
  <> "Returns italicised version of str.\n"
  <> "To be used for AxesLabel etc.\n"
  <> "Automatically threads over lists."
);


Italicised[str_] := Style[str, Italic];


SetAttributes[Italicised, Listable];


(* ::Subsubsection:: *)
(*ModuleSymbol*)


ModuleSymbol::usage = (
  "ModuleSymbol[name, num]\n"
  <> "Builds the symbol with name name and $ModuleNumber num."
);


ModuleSymbol[name_String, num_Integer] :=
  Symbol @ StringJoin[name, "$", ToString[num]];


(* ::Subsubsection:: *)
(*NoExtrapolation*)


(* ::Text:: *)
(*https://mathematica.stackexchange.com/a/30946*)


NoExtrapolation::usage = (
  "NoExtrapolation\n"
  <> "Undocumented option for preventing extrapolation "
  <> "in InterpolatingFunction etc."
);


NoExtrapolation = (
  "ExtrapolationHandler" -> {Indeterminate &, "WarningMessage" -> False}
);


(* ::Subsubsection:: *)
(*PlotOptions*)


PlotOptions::usage = (
  "PlotOptions[type]\n"
  <> "Returns sensible options for a plot of type type, "
  <> "either Axes or Frame."
);


PlotOptions[Axes] = {
  AxesLabel -> Italicised @ {"x", "y"},
  LabelStyle -> Directive[Black, 16]
};


PlotOptions[Frame] = {
  AspectRatio -> Automatic,
  FrameLabel -> Italicised @ {"x", "y"},
  LabelStyle -> Directive[Black, 16],
  RotateLabel -> False
};


(* ::Subsubsection:: *)
(*PreciseOptions*)


PreciseOptions::usage = (
  "Precise[maxStepFrac (def 1/512)]\n"
  <> "Returns options for high precision in NDSolve.\n"
  <> "These may be overridden by options called before PreciseOptions."
);


PreciseOptions[maxStepFrac : _?NumericQ : 1/512] := {
  MaxStepFraction -> maxStepFrac,
  MaxSteps -> Infinity,
  AccuracyGoal -> 16,
  PrecisionGoal -> 16,
  WorkingPrecision -> 2 MachinePrecision
};


(* ::Subsubsection:: *)
(*PrettyString*)


PrettyString::usage = (
  "PrettyString[s1 -> sp1, s2 -> sp2, ...][expr]\n"
  <> "Performs string replacements {s1 -> sp1, s2 -> sp2, ...} "
  <> "on all string subparts of expr.\n"
  <> "To be used to make pretty output with readable input."
);


PrettyString[ruleSeq___Rule][expr_] :=
  expr /. s_String :> StringReplace[s, {ruleSeq}];


(* ::Subsubsection:: *)
(*ReInterpolate*)


ReInterpolate::usage = (
  "ReInterpolate[iFun, num (def 1024)]\n"
  <> "Re-interpolates interpolating function iFun using num uniform subintervals "
  <> "if iFun has more than num subintervals.\n"
  <> "To be used to save space.\n"
  <> "Automatically threads over lists."
);


ReInterpolate[iFun_InterpolatingFunction, num : _?NumericQ : 1024] :=
  (* If more than num subintervals: *)
  If[Length @ iFun["Grid"] - 1 > num,
    (* Re-interpolate *)
    Interpolation @ Table[
      {x, iFun[x]}
    , {x, Subdivide[DomainStart[iFun], DomainEnd[iFun], num]}],
    (* Otherwise return original interpolating function *)
    iFun
  ];


SetAttributes[ReInterpolate, Listable];


(* ::Subsubsection:: *)
(*SeekFirstRootBisection*)


SeekFirstRootBisection::usage = (
  "SeekFirstRootBisection[fun, {xMin, xMax}, num (def 1000), tol (def 10^-6), "
  <> "opts (def N -> True, \"ReturnIterations\" -> False, "
  <> "\"ToleranceType\" -> \"Absolute\")]\n"
  <> "Seeks a numerical root of fun[x] == 0, using the bisection algorithm "
  <> "over the subinterval {x0, x1} corresponding to the first sign-change "
  <> "as determined by num subdivisions of the interval {xMin, xMax}, "
  <> "with tolerance tol over {x0, x1}.\n"
  <> "Default option N -> True specifies that initial endpoints "
  <> "should be numericised (for speed) before performing the iteration.\n"
  <> "Default option \"ReturnIterations\" -> \"False\" "
  <> "specifies that the number of iterations should not be returned.\n"
  <> "Default option \"ToleranceType\" -> \"Absolute\" "
  <> "specifies absolute tolerance. "
  <> "Use \"ToleranceType\" -> \"Relative\" "
  <> "for tolerance relative to the interval {x0, x1}."
);


SeekFirstRootBisection[
  fun_,
  {xMin_?NumericQ, xMax_?NumericQ},
  num : Except[_?OptionQ, _?NumericQ] : 1000,
  tol : Except[_?OptionQ, _?NumericQ] : 10^-6,
  opts : OptionsPattern @ {
    N -> True,
    "ReturnIterations" -> False,
    "ToleranceType" -> "Absolute"
  }
] :=
  Module[{xValues, firstSignChangeIndex, x0, x1},
    (* Determine first sign change subinterval *)
    xValues = (
      Subdivide[xMin, xMax, num]
      // If[TrueQ @ OptionValue[N], N, Identity]
    );
    firstSignChangeIndex = (
      First @ FirstPosition[
        Differences @ Sign[fun /@ xValues],
        signChange_ /; signChange != 0
      ]
    );
    x0 = xValues[[firstSignChangeIndex]];
    x1 = xValues[[firstSignChangeIndex + 1]];
    (* Pass to SeekRootBisection *)
    SeekRootBisection[
      fun,
      {x0, x1},
      tol,
      # -> OptionValue[#] & /@ {
        N,
        "ReturnIterations",
        "ToleranceType"
      }
    ]
  ];


(* ::Subsubsection:: *)
(*SeekRoot*)


SeekRoot::usage = (
  "SeekRoot[fun, {xMin, xMax}, num (def 1000), opts (def N -> True)]\n"
  <> "Seeks a numerical root of fun[x] == 0, using FindRoot, "
  <> "over the search interval {xMin, xMax} "
  <> "with num initial subdivisions to determine the best initial guess.\n"
  <> "Default option N -> True specifies that initial subdivision points "
  <> "should be numericised (for speed) before fun is applied.\n\n"
  <>
  "SeekRoot[fun, xInit]\n"
  <> "Uses initial guess xInit."
);


SeekRoot[
  fun_,
  {xMin_?NumericQ, xMax_?NumericQ},
  num : Except[_?OptionQ, _?NumericQ] : 1000,
  opts : OptionsPattern @ {N -> True}
] :=
  Module[{xInit},
    xInit = First @ TakeSmallestBy[
      Subdivide[xMin, xMax, num] // If[TrueQ @ OptionValue[N], N, Identity],
      Abs @* fun,
      1,
      ExcludedForms -> Except[_?NumericQ] (* exclude infinities *)
    ];
    First @ FindRoot[fun, {xInit}]
  ];


SeekRoot[fun_, xInit_?NumericQ] := First @ FindRoot[fun, {xInit}];


(* ::Subsubsection:: *)
(*SeekRootBisection*)


SeekRootBisection::usage = (
  "SeekRootBisection[fun, {xMin, xMax}, tol (def 10^-6), "
  <> "opts (def N -> True, \"ReturnIterations\" -> False, "
  <> "\"ToleranceType\" -> \"Absolute\")]\n"
  <> "Seeks a numerical root of fun[x] == 0, using the bisection algorithm, "
  <> "over the search interval {xMin, xMax} with tolerance tol.\n"
  <> "Default option N -> True specifies that initial endpoints "
  <> "should be numericised (for speed) before performing the iteration.\n"
  <> "Default option \"ReturnIterations\" -> \"False\" "
  <> "specifies that the number of iterations should not be returned.\n"
  <> "Default option \"ToleranceType\" -> \"Absolute\" "
  <> "specifies absolute tolerance. "
  <> "Use \"ToleranceType\" -> \"Relative\" "
  <> "for tolerance relative to the interval {xMin, xMax}."
);


SeekRootBisection[
  fun_,
  {xMin_?NumericQ, xMax_?NumericQ},
  tol : Except[_?OptionQ, _?NumericQ] : 10^-6,
  opts : OptionsPattern @ {
    N -> True,
    "ReturnIterations" -> False,
    "ToleranceType" -> "Absolute"
  }
] :=
  Module[{absTol, num, xLeft, xRight, xMid, xSol},
    (* Tolerance in x *)
    absTol = tol;
    If[OptionValue["ToleranceType"] =!= "Absolute",
      absTol *= xMax - xMin;
    ];
    (* Initialise *)
    num = 1;
    {xLeft, xRight} = (
      MinMax @ {xMin, xMax} // If[TrueQ @ OptionValue[N], N, Identity]
    );
    While[num < $RecursionLimit,
      (* New midpoint *)
      xMid = Way[xLeft, xRight];
      (* Break if within tolerance *)
      If[(xRight - xLeft) / 2 < absTol,
        xSol = xMid;
        Break[];
      ];
      (* New interval *)
      If[fun[xLeft] fun[xMid] > 0,
        xLeft = xMid,
        xRight = xMid
      ];
      (* Increment iteration count *)
      num += 1;
    ];
    (* Return *)
    If[TrueQ @ Not @ OptionValue["ReturnIterations"],
      xSol,
      {xSol, num}
    ]
  ];


(* ::Subsubsection:: *)
(*SortByPhi*)


SortByPhi::usage = (
  "SortByPhi[list]\n"
  <> "Sorts the list of (planar) points list by azimuthal angle."
);


SortByPhi = SortBy[Mod[ArcTan @@ #, 2 Pi] &]


(* ::Subsubsection:: *)
(*UniformRange*)


UniformRange::usage = (
  "UniformRange[start, end, maxStep]\n"
  <> "Subdivides the interval from start to end uniformly "
  <> "with step at most maxStep."
);


UniformRange[start_?NumericQ, end_?NumericQ, maxStep_?NumericQ] :=
  Subdivide[start, end, (end - start) / maxStep // N // Ceiling]


(* ::Subsubsection:: *)
(*ToName*)


ToName::usage = (
  "ToName[number]\n"
  <> "Converts real number to a string for a file name, "
  <> "with decimal places replaced by underscores."
);


ToName[number_Integer | number_Real] :=
  StringReplace[ToString @ DecimalForm[number], "." -> "_"];


ToName[number_?NumericQ] /; Im[number] == 0 := ToName[Re @ number // N];


(* ::Subsubsection:: *)
(*Way*)


Way::usage = (
  "Way[start, end, prop (def 1/2)]\n"
  <> "Returns proportion prop of the way from start to end.\n"
  <> "The default for prop is 1/2 (halfway)."
);


Way[start_, end_, prop_ : 1/2] := (1 - prop) start + (prop) end;


(* ::Subsubsection:: *)
(*End of private scope*)


End[];


(* ::Subsection:: *)
(*Protect definitions*)


Protect["Conway`*"];


(* ::Subsection:: *)
(*End of package*)


EndPackage[];
