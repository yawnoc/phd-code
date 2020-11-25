(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< Curvilinear`
<< LaplaceYoung`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "general"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Section:: *)
(*Figure: contact angle against a wall (capillary-contact-angle)*)


Module[
  {
    gamma,
    xMin, xMax,
    tStart, tEnd,
    tMin, tMax,
    wallThickness,
    tangentLineLength, angleMarkerRadius,
    normalVectorVerticalPosition, normalVectorLength,
    textStyle, textStyleGreek,
    dummyForTrailingCommas
  },
  (* Constants *)
  gamma = 35 Degree;
  {xMin, xMax} = {0, 1.8};
  (* Endpoints for parameter T *)
  tStart = HHalfPlane[gamma];
  tEnd = SeekRoot[XHalfPlane[gamma][#] - xMax &, {0, tStart}, 10] // Quiet;
  (* Plot range for T (vertical coordinate) *)
  tMin = 0;
  tMax = 1.2 tStart;
  (* Plotting constants *)
  wallThickness = 0.1;
  tangentLineLength = 0.5 (tMax - tMin);
  angleMarkerRadius = 0.4 tangentLineLength;
  normalVectorVerticalPosition = Way[tMin, tMax, 0.45];
  normalVectorLength = 0.4 (tMax - tMin);
  (* Diagram *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStyleGreek = Style[#, LabelSize["LabelOmega"]] & @* LaTeXStyle;
  Show[
    EmptyFrame[{xMin, xMax}, {tStart, tEnd}
      , Frame -> None
      , ImageSize -> 0.6 ImageSizeTextWidth
      , PlotRange -> All
      , PlotRangeClipping -> False
    ],
    (* Contact angle *)
    Graphics @ {Directive[Dashed, GeneralStyle["DefaultThick"]],
      Line @ {
        {xMin, tStart},
        {xMin, tStart} + XYPolar[tangentLineLength, -Pi/2 + gamma]
      }
    },
    Graphics @ {GeneralStyle["DefaultThick"],
      Circle[
        {xMin, tStart},
        angleMarkerRadius,
        {-Pi/2, -Pi/2 + gamma}
      ]
    },
    Graphics @ {
      Text[
        "\[Gamma]" // textStyle
        , {xMin, tStart} + XYPolar[angleMarkerRadius, -Pi/2 + gamma/2]
        , {-0.4, 0.7}
      ]
    },
    (* Liquid *)
    ParametricPlot[
      {XHalfPlane[gamma][t], t}
      , {t, tStart, tEnd}
      , PlotStyle -> Directive[Black, GeneralStyle["Thick"]]
    ],
    Graphics @ {
      Text[
        "liquid" // textStyle
        , {Way[xMin, xMax, 1/5], Way[tMin, tMax, 18/100]}
      ]
    },
    Graphics @ {
      Line @ {{xMin - wallThickness, tMin}, {xMax, tMin}}
    },
    Graphics @ {
      Text[
        "\[CapitalOmega]" // textStyleGreek
        , {Way[xMin - wallThickness/2, xMax], tMin}
        , {0, 0.8}
      ]
    },
    (* Outward normal *)
    Graphics @ {Directive[GeneralStyle["Thick"], Arrowheads[Medium]],
      Arrow @ {
        {xMin - wallThickness, normalVectorVerticalPosition},
        {xMin - wallThickness - normalVectorLength, normalVectorVerticalPosition}
      }
    },
    Graphics @ {
      Text[
        Embolden["n"] // textStyle
        , {xMin - wallThickness - normalVectorLength, normalVectorVerticalPosition}
        , {0, -1.2}
      ]
    },
    (* Wall *)
    Graphics @ {BoundaryTracingStyle["Wall"],
      Rectangle[{xMin - wallThickness, tMin}, {0, tMax}]
    },
    Graphics @ {
      Text[
        "wall" // textStyle
        , {xMin - wallThickness/2, tMax}
        , {0, -1.2}
      ]
    },
    Graphics @ {
      Text[
        "\[PartialD]\[CapitalOmega]" // textStyleGreek
        , {xMin - wallThickness/2, tMin}
        , {0, 0.7}
      ]
    },
    (* Air *)
    Graphics @ {
      Text[
        "air" // textStyle
        , {Way[xMin, xMax, 1/2], Way[tMin, tMax, 65/100]}
      ]
    },
    {}
  ]
] // Ex["capillary-contact-angle.pdf"]
