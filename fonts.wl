(* ::Package:: *)

(* ::Section:: *)
(*List of installed fonts*)


(* ::Text:: *)
(*See https://mathematica.stackexchange.com/a/18945*)
(*Output needs to be tweaked on Mac.*)


ReplaceAll[
  FE`Evaluate @ FEPrivate`GetPopupList["MenuListFonts"],
  (string_ -> family_) :> Style[string, FontFamily -> family]
]
