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


(* RUN ON WINDOWS FOR ITALICISED OMEGA *)
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
      , PlotRangePadding -> {Automatic, Scaled[0.13]}
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
        , If[$OperatingSystem == "Windows", {-0.85, 0.2}, {-0.4, 0.7}]
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
        Italicise["\[CapitalOmega]"] // textStyleGreek
        , {Way[xMin - wallThickness/2, xMax], tMin}
        , If[$OperatingSystem == "Windows", {0, 0.235}, {0.8}]
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
        SeparatedRow[
          If[$OperatingSystem == "Windows",
            StringRepeat["\[NegativeThickSpace]", 9] <> "\[NegativeVeryThinSpace]" // textStyle,
            ""
          ]
        ] @@ {
          "\[PartialD]" // textStyleGreek,
          Italicise["\[CapitalOmega]"] // textStyleGreek
        }
        , {xMin - wallThickness/2, tMin}
        , If[$OperatingSystem == "Windows", {0, 0.9}, {0, 0.7}]
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


(* ::Section:: *)
(*Figure: modification of a wedge (capillary-wedge-modification-*)*)


Module[
  {
    rMaxWall, alpha,
    xMaxWall, yMaxWall,
    rJoin, rCentreLine, rTruncation, xyTruncation,
      truncationPhi, truncationR,
    orthogonalityMarkerLength, thicknessCorrection,
    normalVectorLength, normalVector,
    textStyle,
    options,
    dummyForTrailingCommas
  },
  (* Geometry *)
  rMaxWall = 3;
  alpha = 30 Degree;
  {xMaxWall, yMaxWall} = XYPolar[rMaxWall, alpha];
  (* Truncation boundary *)
  rJoin = 0.4 rMaxWall;
  rCentreLine = 0.3 rMaxWall;
  rTruncation = Interpolation[
    {
      {alpha, rJoin},
      {0, rCentreLine},
      {-alpha, rJoin},
      Nothing
    }
    , InterpolationOrder -> 2
  ];
  xyTruncation[phi_] := XYPolar[rTruncation[phi], phi];
    truncationPhi = 1/3 alpha;
    truncationR = rTruncation[truncationPhi];
  (* Normal vector *)
  orthogonalityMarkerLength = 0.04 rMaxWall;
  thicknessCorrection = 1.07;
  normalVectorLength = 0.3 rMaxWall;
  normalVector[position_, direction_, text_, offset_] := {
    (* Orthogonality marker *)
    Line[orthogonalityMarkerLength {
      {1, 0},
      {1, thicknessCorrection},
      {0, thicknessCorrection}
    }]
      // RotationTransform[ArcTan @@ direction, {0, 0}]
      // TranslationTransform[position]
    ,
    (* Arrow *)
    Directive[GeneralStyle["Thick"], Arrowheads[Medium]],
    Arrow @ {{position, position + normalVectorLength * Normalize[direction]}},
    Text[
      text // textStyle
      , position + normalVectorLength * Normalize[direction]
      , offset
    ]
  };
  (* Text functions *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  (* Options *)
  options = {Nothing
    , ImageSize -> 0.4 * 0.7 ImageSizeTextWidth
    , PlotRange -> {{0, All}, All}
    , PlotRangePadding -> Scaled[0.02]
  };
  (* Make diagrams *)
  {
    (* Original *)
    Show[
      (* Unit normal *)
      Graphics @ {
        normalVector[
          1/4 {xMaxWall, yMaxWall},
          XYPolar[1, Pi/2 + alpha],
          Embolden["n"]
          , {-2.5, -0.5}
        ]
      },
      (* Wedge walls *)
      Graphics @ {BoundaryTracingStyle["Wall"],
        Line @ {
          {xMaxWall, +yMaxWall},
          {0, 0},
          {xMaxWall, -yMaxWall}
        }
      },
      (* Solution *)
      Graphics @ {
        Text[
          Italicise["T"] // textStyle
          , {0.6 xMaxWall, 0}
        ]
      },
      {}
      , options
    ]
      // Ex["capillary-wedge-modification-original.pdf"]
    ,
    (* Modified *)
    Show[
      (* Unit normal *)
      Graphics @ {
        normalVector[
          xyTruncation[truncationPhi],
          Plus[
            rTruncation'[truncationPhi] ARPolar[truncationR, truncationPhi],
            truncationR APhiPolar[truncationR, truncationPhi]
          ] // Cross,
          Subscript[Embolden["n"], "m"]
          , {-1, -1.15}
        ]
      },
      (* Wedge walls *)
      Graphics @ {BoundaryTracingStyle["Wall"],
        Line @ {XYPolar[rJoin, +alpha], {xMaxWall, +yMaxWall}},
        Line @ {XYPolar[rJoin, -alpha], {xMaxWall, -yMaxWall}},
        {}
      },
      ParametricPlot[
        xyTruncation[phi]
        , {phi, alpha, -alpha}
        , PlotPoints -> 2
        , PlotStyle -> BoundaryTracingStyle["Wall"]
      ],
      (* Solution *)
      Graphics @ {
        Text[
          Subscript[Italicise["T"], "\[NegativeThinSpace]\[NegativeVeryThinSpace]m"] // textStyle
          , {Way[rCentreLine, xMaxWall, 0.55], 0}
        ]
      },
      {}
      , options
    ]
      // Ex["capillary-wedge-modification-modified.pdf"]
  }
]


(* ::Section:: *)
(*Figure: serrated approximation geometry for a modified boundary (capillary-serrated-approximation-geometry-modified)*)


Module[
  {
    theta,
    tracedLongLength, tracedShortLength,
    tracedLongDisplacement, tracedShortDisplacement,
    tracedXPeriod,
    origin,
    contourAnchor,
    thetaM,
    angleMarkerRadiusThetaLeft,
    angleMarkerRadiusThetaRight,
    angleMarkerRadiusThetaM,
    gradientVectorLength, normalVectorLength,
    arrowStyle,
    orthogonalityMarkerLength, orthogonalityMarker,
    textStyleStraight, textStyleLabel, textStyleLabelBigger,
      normalAnchor, normalVector, normalTip,
      normalMAnchor, normalMVector, normalMTip,
      thetaAnchor,
      thetaMAnchor,
      modifiedAnchorLeft, modifiedAnchorRight,
      modifiedStart, modifiedEnd,
      contourStart, contourEnd,
      gradientAnchor, gradientTip,
    dummyForTrailingCommas
  },
  (* Geometry *)
  theta = 50 Degree;
  {tracedLongLength, tracedShortLength} = {2, 0.9};
  tracedLongDisplacement = XYPolar[tracedLongLength, theta];
  tracedShortDisplacement = XYPolar[tracedShortLength, -theta];
  tracedXPeriod = tracedLongDisplacement + tracedShortDisplacement // First;
  origin = {0, 0};
  contourAnchor = origin + tracedLongDisplacement + tracedShortDisplacement;
  thetaM = ArcTan @@ (tracedLongDisplacement + tracedShortDisplacement);
  (* Etc. *)
  angleMarkerRadiusThetaLeft = 0.23 tracedXPeriod;
  angleMarkerRadiusThetaRight = 0.12 tracedXPeriod;
  angleMarkerRadiusThetaM = 0.23 tracedXPeriod;
  gradientVectorLength = 0.45 tracedXPeriod;
  normalVectorLength = 0.4 tracedXPeriod;
  arrowStyle = Directive[GeneralStyle["Thick"], Arrowheads[0.04]];
  orthogonalityMarkerLength = 0.13 normalVectorLength;
  orthogonalityMarker[position_, angle_] := {
    Line[orthogonalityMarkerLength {{1, 0}, {1, 1}, {0, 1}}]
      // RotationTransform[angle, {0, 0}]
      // TranslationTransform[position]
  };
  (* Text styles *)
  textStyleStraight = Style[#, LabelSize["Straight"]] & @* LaTeXStyle;
  textStyleLabel = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStyleLabelBigger = Style[#, LabelSize["Label"] + 1] & @* LaTeXStyle;
  (* Make diagram *)
  Show[
    thetaAnchor = contourAnchor - 2 {First[tracedShortDisplacement], 0};
    (* Angle between traced boundaries and contour *)
    Graphics @ {
      Circle[thetaAnchor, angleMarkerRadiusThetaLeft, {0, theta}],
      Circle[contourAnchor, angleMarkerRadiusThetaRight, {Pi, Pi - theta}],
      {}
    },
    Graphics @ {
      Text[
        "\[Theta]" // textStyleLabelBigger (* optical illusion *)
        , thetaAnchor + XYPolar[angleMarkerRadiusThetaLeft, Way[thetaM, theta]]
        , {-2, -0.2}
      ]
    },
    Graphics @ {
      Text[
        "\[Theta]" // textStyleLabelBigger (* optical illusion *)
        , contourAnchor + XYPolar[angleMarkerRadiusThetaRight, Pi - theta/2]
        , {1.4, -0.2}
      ]
    },
    (* Serrated path normal vector *)
    normalAnchor = origin + 0.85 tracedLongDisplacement;
    normalVector = XYPolar[normalVectorLength, theta + Pi/2];
    normalTip = normalAnchor + normalVector;
    Graphics @ {
      orthogonalityMarker[normalAnchor, theta]
    },
    Graphics @ {arrowStyle,
      Arrow @ {normalAnchor, normalTip}
    },
    Graphics @ {
      Text[
        Embolden["n"] // textStyleLabel
        , normalTip
        , {1, -0.9}
      ]
    },
    (* Modified boundary normal vector *)
    normalMAnchor = (
      contourAnchor
      - 1/2 tracedShortDisplacement
      + 0.3 tracedXPeriod XYPolar[1, thetaM]
    );
    normalMVector = XYPolar[normalVectorLength, thetaM + Pi/2];
    normalMTip = normalMAnchor + normalMVector;
    Graphics @ {
      orthogonalityMarker[normalMAnchor, thetaM]
    },
    Graphics @ {arrowStyle,
      Arrow @ {normalMAnchor, normalMTip}
    },
    Graphics @ {
      Text[
        Subscript[Embolden["n"], "m"] // textStyleLabel
        , normalMTip
        , {-0.35, -0.9}
      ]
    },
    (* Angle between modified boundary and contour *)
    thetaMAnchor = (
      contourAnchor
      - 1/2 tracedShortDisplacement
      + {Cot[thetaM], 1} 1/2 Last[tracedShortDisplacement]
    );
    Graphics @ {
      Circle[thetaMAnchor, angleMarkerRadiusThetaM, {0, thetaM}]
    },
    Graphics @ {
      Text[
        Subscript["\[Theta]", "\[NegativeVeryThinSpace]\[NegativeVeryThinSpace]m"] // textStyleLabel
        , thetaMAnchor + XYPolar[angleMarkerRadiusThetaM, thetaM/2]
        , {-1.4, -0.25}
      ]
    },
    (* Modified boundary *)
    modifiedAnchorLeft = origin - tracedShortDisplacement / 2;
    modifiedAnchorRight = origin + 2 tracedLongDisplacement + 3/2 tracedShortDisplacement;
    modifiedStart = Way[modifiedAnchorLeft, modifiedAnchorRight, -0.15];
    modifiedEnd = Way[modifiedAnchorLeft, modifiedAnchorRight, +1.15];
    Graphics @ {BoundaryTracingStyle["Wall"],
      Line @ {modifiedStart, modifiedEnd}
    },
    Graphics @ {Darker @ Darker[Gray],
      Text[
        "modified boundary" // textStyleStraight
        , modifiedEnd
        , {0.95, 0.8}
        , XYPolar[1, thetaM]
      ]
    },
    (* Serrated path *)
    Graphics @ {BoundaryTracingStyle["Traced"],
      Line @ {
        origin - tracedShortDisplacement,
        origin,
        origin + tracedLongDisplacement,
        origin + tracedLongDisplacement + tracedShortDisplacement,
        origin + 2 tracedLongDisplacement + tracedShortDisplacement,
        origin + 2 tracedLongDisplacement + 2 tracedShortDisplacement,
        Nothing
      }
    },
    Graphics @ {
      Text[
        "serrated path" // textStyleStraight
        , origin
        , {-1.27, -0.6}
      ]
    },
    (* T-contour *)
    contourStart = contourAnchor + tracedXPeriod {-1.4, 0};
    contourEnd = contourAnchor + tracedXPeriod {1.2, 0};
    Graphics @ {BoundaryTracingStyle["Contour"],
      Line @ {contourStart, contourEnd}
    },
    Graphics @ {
      Text[
        Row @ {Italicise["T"], "\[Hyphen]contour"} // textStyleStraight
        , contourEnd
        , {1.4, -1.1}
      ]
    },
    (* Gradient vector *)
    gradientAnchor = thetaAnchor - 0.5 {tracedXPeriod, 0};
    gradientTip = gradientAnchor + {0, gradientVectorLength};
    Graphics @ {
      orthogonalityMarker[gradientAnchor, 0]
    },
    Graphics @ {arrowStyle,
      Arrow @ {gradientAnchor, gradientTip}
    },
    Graphics @ {
      (* Fake Del *)
      {
        RegularPolygon[3]
          // Scale[#, 0.095 * 0.77] &
          // Rotate[#, Pi] &
          // Translate[#, gradientVectorLength {-0.336, 1.843}] &
      },
      {White,
        RegularPolygon[3]
          // Scale[#, 0.045 * 0.77] &
          // Rotate[#, Pi] &
          // Translate[#, gradientVectorLength {-0.329, 1.847}] &
      },
      Text[
        Italicise["T"] // textStyleLabel
        , gradientTip
        , {-1, -1}
      ],
      {}
    },
    {}
    , ImageSize -> 0.9 ImageSizeTextWidth
  ]
] // Ex["capillary-serrated-approximation-geometry-modified.pdf"]
