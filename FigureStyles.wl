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
  BoundaryTracingStyle,
  SquigglyArrow,
  {}
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
  "VeryThick" -> AbsoluteThickness[5],
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
  "Wall" -> Directive[
    GeneralStyle["VeryThick"],
    Gray
  ],
  Automatic -> Black
][type];


BoundaryTracingStyle[typeSeq__String] :=
  Directive @@ BoundaryTracingStyle /@ {typeSeq};


(* ::Subsubsection:: *)
(*SquigglyArrow*)


SquigglyArrow[{xBase_, yBase_}, phi_: 0, size_: 1] :=
  Module[{yList, listLength, xList, pointList},
    (* List of y-coordinates (for x-coordinates 0, ..., 1) *)
    (*
      NOTE:
      - Nonzero y-coordinate in the ConstantArrays is an optical correction
        to make the straight parts of the arrow appear in line.
      - Greater number of elements in the second ConstantArray is a correction
        to make the arrow appear balanced when the arrowhead is included.
    *)
    yList =
      {
        ConstantArray[-0.1, 5],
        1, 0, -1, 0, 1, 0, -1,
        ConstantArray[+0.1, 7],
        {}
      }
        // 0.15 # &
        // Flatten;
    (* List of x-coordinates *)
    listLength = Length[yList];
    xList = Subdivide[0, 1, listLength - 1];
    (* Make B-spline arrow *)
    pointList = {xList, yList} // Transpose;
    Arrow @ BSplineCurve[
      pointList
        // ScalingTransform @ {size, size}
        // RotationTransform[phi, {0, 0}]
        // TranslationTransform @ {xBase, yBase}
    ]
  ];


(* ::Subsubsection:: *)
(*End of private scope*)


End[];


(* ::Subsection:: *)
(*Protect definitions*)


Protect["FigureStyles`*"];


(* ::Subsection:: *)
(*End of package*)


EndPackage[];
