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
    Graphics @ {
      Line @ {
        {xMin, tStart},
        {xMin, tStart} + XYPolar[tangentLineLength, -Pi/2 + gamma]
      }
    },
    Graphics @ {
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


(* ::Section:: *)
(*Figure: serrated approximation geometry (capillary-serrated-approximation-geometry)*)


Module[
  {
    theta,
    edgeHalfWidth, edgeHalfHeight,
    orthogonalityMarkerLength, angleMarkerLength,
    normalVectorLength, gradientVectorLength,
    arrowStyle,
    textStyleStraight, textStyleLabel,
    dummyForTrailingCommas
  },
  (* Geometry *)
  theta = 30 Degree;
  edgeHalfWidth = 1;
  edgeHalfHeight = edgeHalfWidth Tan[theta];
  (* Etc. *)
  orthogonalityMarkerLength = 0.15 normalVectorLength;
  angleMarkerLength = 0.35 edgeHalfWidth;
  normalVectorLength = 0.9 edgeHalfWidth;
  gradientVectorLength = 1.2 edgeHalfWidth;
  arrowStyle = Directive[GeneralStyle["Thick"], Arrowheads[0.04]];
  (* Text styles *)
  textStyleStraight = Style[#, LabelSize["Straight"]] & @* LaTeXStyle;
  textStyleLabel = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  (* Make diagram *)
  Show[
    (* Serrated path normal vectors *)
    Graphics @ {
      Line[orthogonalityMarkerLength {{0, 1}, {1, 1}, {1, 0}}]
        // Rotate[#, theta, {0, 0}] &
        // Translate[#, {-edgeHalfWidth, 0}] &,
      Line[orthogonalityMarkerLength {{0, 1}, {-1, 1}, {-1, 0}}]
        // Rotate[#, -theta, {0, 0}] &
        // Translate[#, {+edgeHalfWidth, 0}] &,
      {}
    },
    Graphics @ {arrowStyle,
      Arrow @ {
        {-edgeHalfWidth, 0},
        {-edgeHalfWidth, 0} + XYPolar[normalVectorLength, Pi/2 + theta]
      },
      Arrow @ {
        {+edgeHalfWidth, 0},
        {+edgeHalfWidth, 0} + XYPolar[normalVectorLength, Pi/2 - theta]
      },
      {}
    },
    Graphics @ {
      Text[
        Embolden["n"] // textStyleLabel
        , {-edgeHalfWidth, 0} + XYPolar[normalVectorLength, Pi/2 + theta]
        , {1, -1}
      ],
      Text[
        Embolden["n"] // textStyleLabel
        , {+edgeHalfWidth, 0} + XYPolar[normalVectorLength, Pi/2 - theta]
        , {-1, -1}
      ],
      {}
    },
    (* Angle between normal vector and gradient vector *)
    Graphics @ {
      Circle[{-edgeHalfWidth, 0}, angleMarkerLength, {0, theta}],
      Circle[{+edgeHalfWidth, 0}, angleMarkerLength, {Pi, Pi - theta}],
      {}
    },
    Graphics @ {
      Text[
        "\[Theta]" // textStyleLabel
        , {-edgeHalfWidth, 0} + XYPolar[angleMarkerLength, theta/2]
        , {-2.2, -0.1}
      ],
      Text[
        "\[Theta]" // textStyleLabel
        , {+edgeHalfWidth, 0} + XYPolar[angleMarkerLength, Pi - theta/2]
        , {1.7, -0.1}
      ],
      {}
    },
    (* Serrated path *)
    Graphics @ {BoundaryTracingStyle["Traced"],
      Line @ {
        {-3 edgeHalfWidth, 0},
        {-2 edgeHalfWidth, -edgeHalfHeight},
        {0, edgeHalfHeight},
        {2 edgeHalfWidth, -edgeHalfHeight},
        Way[
          {2 edgeHalfWidth, -edgeHalfHeight},
          {3 edgeHalfWidth, 0}
          , 0.7
        ],
        Nothing
      }
    },
    Graphics @ {
      Text[
        "serrated path" // textStyleStraight
        , {2 edgeHalfWidth, -edgeHalfHeight}
        , {0, 1}
      ]
    },
    (* T-contour *)
    Graphics @ {BoundaryTracingStyle["Contour"],
      Line @ {
        {-3.3 edgeHalfWidth, 0},
        {3.1 edgeHalfWidth, 0}
      }
    },
    Graphics @ {
      Text[
        Row @ {Italicise["T"], "\[Hyphen]contour"} // textStyleStraight
        , {2.5 edgeHalfWidth, 0}
        , {0, -1.2}
      ]
    },
    (* Gradient vector *)
    Graphics @ {
      Line[orthogonalityMarkerLength {{0, 1}, {-1, 1}, {-1, 0}}]
    },
    Graphics @ {arrowStyle,
      Arrow @ {{0, 0}, {0, gradientVectorLength}}
    },
    Graphics @ {
      (* Fake Del *)
      {
        RegularPolygon[3]
          // Scale[#, 0.095] &
          // Rotate[#, Pi] &
          // Translate[#, gradientVectorLength {-0.06, 0.911}] &
      },
      {White,
        RegularPolygon[3]
          // Scale[#, 0.045] &
          // Rotate[#, Pi] &
          // Translate[#, gradientVectorLength {-0.053, 0.915}] &
      },
      Text[
        Italicise["T"] // textStyleLabel
        , {0, gradientVectorLength}
        , {-1, -1}
      ],
      {}
    },
    {}
    , ImageSize -> 0.9 ImageSizeTextWidth
  ]
] // Ex["capillary-serrated-approximation-geometry.pdf"]


(* ::Section:: *)
(*Figure: wall with grooves (capillary-wall-with-grooves)*)


Module[
  {
    sigma, lambda, phi,
    yMaxGroove, yArrow,
    angleMarkerRadius,
    arrowStyle,
    textStyleLabel,
    dummyForTrailingCommas
  },
  (* Constants *)
  sigma = 2;
  lambda = 4;
  phi = 50 Degree;
  (* More constants *)
  yMaxGroove = sigma/2 Tan[phi];
  yArrow = -0.15 sigma;
  (* Etc. *)
  angleMarkerRadius = 0.22 sigma;
  arrowStyle = Arrowheads @ {-Small, Small};
  textStyleLabel = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  (* Make diagram *)
  Show[
    (* Contour *)
    Graphics @ {
      BoundaryTracingStyle["Contour"],
      GeneralStyle["DefaultThick"],
      Table[
        Line @ {{xGrooveCorner, 0}, {xGrooveCorner + sigma, 0}}
        , {xGrooveCorner, {0, sigma + lambda}}
      ]
    },
    (* Groove spacing *)
    Graphics @ {
      arrowStyle,
      Arrow @ {{sigma, yArrow}, {sigma + lambda, yArrow}},
      Text[
        "\[Lambda]" // textStyleLabel
        , {sigma + lambda/2, yArrow}
        , {0, 1}
      ],
      {}
    },
    (* Groove width *)
    Table[
      Graphics @ {
        arrowStyle,
        Arrow @ {{xGrooveCorner, yArrow}, {xGrooveCorner + sigma, yArrow}},
        Text[
          "\[Sigma]" // textStyleLabel
          , {xGrooveCorner + sigma/2, yArrow}
          , {0, 1}
        ],
        {}
      }
      , {xGrooveCorner, {0, sigma + lambda}}
    ],
    (* Groove angle *)
    Table[
      Graphics @ {
        Circle[{xGrooveCorner, 0}, angleMarkerRadius, {0, phi}],
        Text[
          "\[CurlyPhi]" // textStyleLabel
          , {xGrooveCorner, 0} + XYPolar[angleMarkerRadius, phi/2]
          , {-1.65, -0.5}
        ],
        {}
      }
      , {xGrooveCorner, {0, sigma + lambda}}
    ],
    (* Indented wall *)
    Graphics @ {
      Directive[
        FaceForm[BoundaryTracingStyle["Wall"] // Last],
        EdgeForm[GeneralStyle["DefaultThick"]]
      ],
      Polygon @ {
        {0, 0},
        {sigma/2, yMaxGroove},
        {sigma, 0},
        {sigma + lambda, 0},
        {sigma + lambda + sigma/2, yMaxGroove},
        {2 sigma + lambda, 0},
        {2 sigma + 2 lambda, 0},
        {2 sigma + 2 lambda, 2 yMaxGroove},
        {-lambda, 2 yMaxGroove},
        {-lambda, 0}
      }
    },
    {}
    , ImageSize -> 0.6 ImageSizeTextWidth
    , PlotRange -> {{-0.3 sigma, 2.8 sigma + lambda}, {All, 1.7 yMaxGroove}}
  ]
] // Ex["capillary-wall-with-grooves.pdf"]
