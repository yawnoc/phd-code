(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< FigureStyles`
SetDirectory @ FileNameJoin @ {NotebookDirectory[], "general"}


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Section:: *)
(*Figure: radiation--conduction BVP (radiation-conduction-bvp)*)


Module[
  {
    conductionA, conductionB, conductionPhi,
    conductionXY, conductionNormalPhi,
    bathX, bathY, bathA, bathB,
    textStyle, textStyleGreek,
    radiationArrowClearanceFactor,
    dummyForTrailingCommas
  },
  (* Conduction ellipse *)
  {conductionA, conductionB} = {8, 5};
  conductionPhi = 20 Degree;
  (* Geometry *)
  (* (see <https://math.stackexchange.com/a/990013> for angle of normal) *)
  conductionXY[ang_] := {conductionA Cos[ang], conductionB Sin[ang]};
  conductionNormalPhi[ang_] := ArcTan[conductionB Cos[ang], conductionA Sin[ang]];
  (* Heat bath ellipse *)
  {bathX, bathY} = {-3, -2};
  {bathA, bathB} = {2.5, 1.5};
  (* Diagram *)
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  textStyleGreek = Style[#, LabelSize["LabelOmega"]] & @* LaTeXStyle;
  Show[
    (* Conduction ellipse *)
    Graphics @ {BoundaryTracingStyle["Traced"],
      Circle[{0, 0}, {conductionA, conductionB}]
        // Rotate[#, conductionPhi] &
    },
    Graphics @ {
      Text[
        SeparatedRow["\[ThinSpace]" // textStyle] @@ {
          "conduction" // textStyle,
          "\[CapitalOmega]" // textStyleGreek
        }
        , -{0.3 bathX, 0.8bathY}
      ]
    },
    (* Radiation arrows *)
    radiationArrowClearanceFactor = 1.1;
    Graphics @ {Arrowheads[0.02],
      Table[
        SquigglyArrow[
          radiationArrowClearanceFactor * conductionXY[ang]
          , conductionNormalPhi[ang]
          , 2
        ]
        , {ang, Subdivide[0, 2 Pi, 8]}
      ]
        // Rotate[#, conductionPhi] &
    },
    Graphics @ {
      Text[
        Column[
          {
            "radiation" // textStyle,
            "\[PartialD]\[CapitalOmega]" // textStyleGreek
          }
          , Alignment -> Right
          , Spacings -> 0.3
        ]
        , conductionXY[5/8 Pi] // RotationTransform[conductionPhi]
        , {1.1, -1.15}
      ]
    },
    (* Heat bath ellipse *)
    Graphics @ {Directive[FaceForm @ GrayLevel[0.9], EdgeForm[Black]],
      Disk[{bathX, bathY}, {bathA, bathB}]
    },
    Graphics @ {
      Text[
        "heat" // textStyle
        , {bathX, bathY}
        , {0, -0.1}
      ]
    },
    {}
    , ImageSize -> 0.45 ImageSizeTextWidth
  ]
] // Ex["radiation-conduction-bvp.pdf"]


(* ::Section:: *)
(*Figure: self viewing radiation elements (self-viewing-radiation-elements)*)


Module[
  {
    r, n, rStar, nStar, dStar, theta, thetaStar,
    normalVectorLength, elementWidth,
    angleMarkerRadius, angleMarkerRadiusStar,
    nTip, nStarTip,
    zVector, elementGeneric, element, elementStar,
    sweep,
    centre, eye, canvasNormal, canvasPoint,
    project,
    elementStyle, positionStyle,
    normalVectorStyle, displacementVectorStyle, angleMarkerStyle,
    starred, textStyle,
    dummyForTrailingCommas
  },
  (* Coordinates *)
  r = {-1, 0, 0};
  n = {0, 0, 1} // Normalize;
  rStar = {0.7, 0, 0.3};
  nStar = {-0.4, -0.2, 0.7} // Normalize;
  dStar = r - rStar;
  theta = ArcCos[n . (-dStar) / Norm[dStar]];
  thetaStar = ArcCos[nStar . dStar / Norm[dStar]];
  (* Plotting constants *)
  normalVectorLength = 0.8;
  elementWidth = 0.35;
  angleMarkerRadius = 0.2 normalVectorLength;
  angleMarkerRadiusStar = 0.25 normalVectorLength;
  nTip = r + normalVectorLength * n;
  nStarTip = rStar + normalVectorLength * nStar;
  (* Element coordinates *)
  zVector = {0, 0, 1};
  elementGeneric = elementWidth/2 {{-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}};
  element = (elementGeneric
    // RotationTransform[{zVector, n}]
    // TranslationTransform[r]
  );
  elementStar = (elementGeneric
    // RotationTransform[{zVector, nStar}]
    // TranslationTransform[rStar]
  );
  (* Sweeping transformation (for angle markers) *)
  sweep[origin_, angle_][startVector_, endVector_][vector_] := (
    vector
      // RotationTransform[angle, {startVector, endVector}, {0, 0, 0}]
      // TranslationTransform[origin]
  );
  (* Perspective geometry *)
  centre = Way[r, rStar];
  canvasNormal = {0.15, -1, 0.4};
  canvasPoint = centre + canvasNormal;
  eye = centre + 10 canvasNormal;
  project[point_] :=
    OnePointPerspective[2][eye, canvasNormal, canvasPoint, zVector][point];
  (* Plot styles *)
  elementStyle = Directive[FaceForm[None], EdgeForm[Black]];
  positionStyle = PointSize[Medium];
  normalVectorStyle = Directive[GeneralStyle["Thick"], Arrowheads[0.06]];
  displacementVectorStyle = Arrowheads @ {{0.05, 0.55}};
  angleMarkerStyle = Directive[Black, Thickness[Small]];
  (* Text functions *)
  starred[expr_] := Superscript[expr, "\[NegativeVeryThinSpace]\[FivePointedStar]"];
  textStyle = Style[#, LabelSize["Label"]] & @* LaTeXStyle;
  (* Make scene *)
  Show[
    (* Element dA *)
    Graphics @ {
      positionStyle, Point[project /@ {r}],
        Text[Embolden["r"] // textStyle
          , project[r]
          , {2.7, -0.23}
        ],
      elementStyle, Polygon[project /@ element],
        Text[SeparatedRow[""] @ {"d", Italicise["A"]} // textStyle
          , project[r]
          , {0.5, 1.6}
        ],
      normalVectorStyle, Arrow[project /@ {r, nTip}],
        Text[Embolden["n"] // textStyle
          , project[nTip]
          , {-2.5, 0}
        ],
      {}
    },
    (* Angle theta *)
    ParametricPlot[
      sweep[r, ang][n, -dStar][angleMarkerRadius * n]
        // project
        // Evaluate
      , {ang, 0, theta}
      , PlotPoints -> 2
      , PlotStyle -> angleMarkerStyle
    ],
    Graphics @ {
      Text["\[Theta]" // textStyle
        , sweep[r, theta/2][n, -dStar][angleMarkerRadius * n] // project
        , {-1.5, -0.6}
      ]
    },
    (* Element dA-star *)
    Graphics @ {
      positionStyle, Point[project /@ {rStar}],
        Text[
          Embolden["r"] // starred // textStyle
          , project[rStar]
          , {-1.6, -0.7}
        ],
      elementStyle, Polygon[project /@ elementStar],
        Text[SeparatedRow[""] @ {"d", Italicise["A"] // starred} // textStyle
          , project[rStar]
          , {-0.5, 1.8}
        ],
      normalVectorStyle, Arrow[project /@ {rStar, nStarTip}],
        Text[Embolden["n"] // starred // textStyle
          , project[nStarTip]
          , {-1, -1}
        ],
      {}
    },
    (* Angle theta-star *)
    ParametricPlot[
      sweep[rStar, ang][nStar, dStar][angleMarkerRadiusStar * nStar]
        // project
        // Evaluate
      , {ang, 0, thetaStar}
      , PlotPoints -> 2
      , PlotStyle -> angleMarkerStyle
    ],
    Graphics @ {
      Text["\[Theta]" // starred // textStyle
        , sweep[rStar, theta/2][nStar, dStar][angleMarkerRadiusStar * nStar] // project
        , {0.9, -0.3}
      ]
    },
    (* Displacement d-star *)
    Graphics @ {
      displacementVectorStyle, Arrow[project /@ {rStar, r}],
        Text[
          Embolden["d"] // starred // textStyle
          , project @ Way[rStar, r]
          , {-0.3, 1.3}
        ],
      {}
    },
    {}
    , ImageSize -> 0.5 ImageSizeTextWidth
  ]
] // Ex["self-viewing-radiation-elements.pdf"]
