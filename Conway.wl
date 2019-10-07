(* ::Package:: *)

(* ::Text:: *)
(*Conway.wl (for Wolfram Mathematica 12.0)*)


(* ::Text:: *)
(*A Wolfram Language Package (.wl) with handy functions.*)
(*(C) 2019 Conway*)
(*ABSOLUTELY NO WARRANTY, i.e. "GOD SAVE YOU"*)


(* ::Section:: *)
(*Improve workflow*)


(* ::Text:: *)
(*Disable annoying floating elements (which pop up when I don't need them)*)
(*and make alternative input auto replacements*)
(*(to avoid unreadable \[...] representations in text editors).*)


SetOptions[
  $FrontEndSession,
  CodeAssistOptions -> {
    "FloatingElementEnable" -> False
  },
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
  DomainEnd,
  DomainStart,
  EmptyAxes,
  EmptyFrame,
  Ex,
  Italicised,
  NoExtrapolation,
  PlotOptions,
  PrettyString,
  SeekRoot,
  SeekRootBisection,
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
(*SeekRoot*)


SeekRoot::usage = (
  "SeekRoot[fun, {xMin, xMax}, num (def 1000), opts (def N -> True)]\n"
  <> "Seeks a numerical root of fun[x] == 0, using FindRoot, "
  <> "over the search interval {xMin, xMax} "
  <> "with num initial subdivisions to determine the best initial guess.\n"
  <> "Default option N -> True specifies that initial subdivision points "
  <> "should be numericised (for speed) before fun is applied."
);


SeekRoot[
  fun_,
  {xMin_?NumericQ, xMax_?NumericQ},
  num : Except[_?OptionQ, _?NumericQ] : 1000,
  opts : OptionsPattern @ {N -> True}
] :=
  Module[{xInit, x},
    xInit = First @ TakeSmallestBy[
      Subdivide[xMin, xMax, num] // If[TrueQ @ OptionValue[N], N, Identity],
      Abs @* fun,
      1,
      ExcludedForms -> Except[_?NumericQ] (* exclude infinities *)
    ];
    x /. FindRoot[fun[x], {x, xInit}]
  ];


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
