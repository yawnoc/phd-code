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
  "Point" -> PointSize[Large],
  "Translucent" -> Opacity[0.7],
  Automatic -> Automatic
][type];


(* ::Subsubsection:: *)
(*BoundaryTracingStyle*)


BoundaryTracingStyle[type : (_String | Automatic) : Automatic] := Association[
  "Background" -> GrayLevel[0.7],
  "Contour" -> Directive[AbsoluteDashing @ {1, 4}, Black],
  "ContourImportant" -> Gray,
  "NonViable" -> Directive[GeneralStyle["Translucent"], LightGray],
  "Terminal" -> Directive[AbsoluteDashing @ {12, 6}, Black],
  "Traced" -> Directive[AbsoluteThickness[2.5], Black],
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
