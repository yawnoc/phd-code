(* ::Package:: *)

(* ::Section:: *)
(*Initialisation section (always run this first)*)


(* ::Subsection:: *)
(*Load packages and set export directory*)


SetDirectory @ ParentDirectory @ NotebookDirectory[];
<< Conway`
<< FigureStyles`
SetDirectory @ NotebookDirectory[]


(* ::Subsection:: *)
(*Clean slate*)


ClearAll["Global`*"];


(* ::Section:: *)
(*Figure: radiation BVP (radiation-bvp)*)


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
  textStyle = Style[#, 24] & @* LaTeXStyle;
  textStyleGreek = Style[#, 30] & @* LaTeXStyle;
  Show[
    (* Conduction ellipse *)
    Graphics @ {BoundaryTracingStyle["Traced"],
      Circle[{0, 0}, {conductionA, conductionB}]
        // Rotate[#, conductionPhi] &
    },
    Graphics @ {
      Text[
        SeparatedRow[Invisible["i" // textStyle]] @@ {
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
          , Spacings -> 0.8
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
    , ImageSize -> 360
  ]
] // Ex["radiation-bvp.pdf"]
