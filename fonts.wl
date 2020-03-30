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


(* ::Section:: *)
(*LaTeX fonts*)


(* ::Text:: *)
(*See https://tex.stackexchange.com/q/524996*)
(*and http://mirrors.ctan.org/macros/latex/contrib/unicode-math/unimath-symbols.pdf*)
(*Special unicode symbols need to be used.*)


(* ::Subsection:: *)
(*Symbols defined by unicode-math*)


Module[{fontSize, styleMath, styleRomanItalic},
  fontSize = 36;
  styleMath[string_] := Style[string, fontSize,
    FontFamily -> "Latin Modern Math"
  ];
  styleRomanItalic[string_] := Style[string, fontSize,
    Italic,
    FontFamily -> "Latin Modern Roman"
  ];
  SetAttributes[styleMath, Listable];
  SetAttributes[styleRomanItalic, Listable];
  styleMath @ {
    (* Digits *)
    CharacterRange["0", "9"],
    (* Latin upright *)
    CharacterRange["A", "Z"],
    CharacterRange["a", "z"],
    (* Latin italic *)
    CharacterRange[16^^1D434, 16^^1D44D],
    CharacterRange[16^^1D44E, 16^^1D467],
    (* Latin italic with missing "h" filled in *)
    (*
      NOTE: https://tex.stackexchange.com/questions/524996/#comment1328261_525087
        For the missing "h", blame the Unicode Consortium;
        they already defined the math italic "h" as the Planck constant U+210E
        and didn't fill the hole in the Mathematical Italic block,
        which is a silly decision.
        \[Dash] egreg
     *)
    styleRomanItalic @ CharacterRange["a", "z"],
    (* Greek upright upper and italic lower *)
    CharacterRange[16^^00391, 16^^003A9],
    CharacterRange[16^^1D6FC, 16^^1D71B],
    (* Bold latin upright *)
    CharacterRange[16^^1D400, 16^^1D419],
    CharacterRange[16^^1D41A, 16^^1D433],
    {}
  } // TableForm[#, TableSpacing -> {4, 1}] &
]
