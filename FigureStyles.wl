(* ::Package:: *)

(* ::Text:: *)
(*FigureStyles.wl (for Wolfram Mathematica)*)


(* ::Text:: *)
(*Semantic styling for figures.*)
(*ABSOLUTELY NO WARRANTY, i.e. "GOD SAVE YOU"*)


(* ::Section:: *)
(*Package definitions*)


(* ::Subsection:: *)
(*Start of package*)


BeginPackage["FigureStyles`"];


(* ::Subsection:: *)
(*Clear existing definitions if any*)


Unprotect["FigureStyles`*"];
ClearAll["FigureStyles`*"];
ClearAll["FigureStyles`*`*"];


(* ::Subsection:: *)
(*Mention non-private symbols*)


{
  GeneralStyle,
  BoundaryTracingStyle
};


(* ::Subsection:: *)
(*Private scope*)


(* ::Subsubsection:: *)
(*Start of private scope*)


Begin["`Private`"];


(* ::Subsubsection:: *)
(*GeneralStyle*)


GeneralStyle[type_String : Automatic] := Association[
  "Dashed" -> AbsoluteDashing @ {12, 6},
  "DefaultThick" -> AbsoluteThickness[1.5],
  "Dotted" -> AbsoluteDashing @ {1, 4},
  "Point" -> PointSize[Large],
  "Thick" -> AbsoluteThickness[2.5],
  "Translucent" -> Opacity[0.7],
  Automatic -> Automatic
][type];


(* ::Subsubsection:: *)
(*BoundaryTracingStyle*)


BoundaryTracingStyle[type : (_String | Automatic) : Automatic] := Association[
  "Background" -> Directive[
    GeneralStyle["DefaultThick"],
    GrayLevel[0.75]
  ],
  "Contour" -> Directive[
    GeneralStyle["DefaultThick"],
    GeneralStyle["Dotted"],
    Black
  ],
  "ContourImportant" -> Directive[
    GeneralStyle["Thick"],
    Gray
  ],
  "NonViable" -> Directive[
    GeneralStyle["Translucent"],
    LightGray
  ],
  "Terminal" -> Directive[
    GeneralStyle["DefaultThick"],
    GeneralStyle["Dashed"],
    Black
  ],
  "Traced" -> Directive[
    GeneralStyle["Thick"],
    Black
  ],
  "Unphysical" -> Black,
  "Viable" -> White,
  Automatic -> Black
][type];


BoundaryTracingStyle[typeSeq__String] :=
  Directive @@ BoundaryTracingStyle /@ {typeSeq};


(* ::Subsubsection:: *)
(*End of private scope*)


End[];


(* ::Subsection:: *)
(*Protect definitions*)


Protect["FigureStyles`*"];


(* ::Subsection:: *)
(*End of package*)


EndPackage[];
