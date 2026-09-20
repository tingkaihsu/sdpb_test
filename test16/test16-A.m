(* ::Package:: *)

sourceDirectory = If[
  StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
  DirectoryName[ExpandFileName[$InputFileName]],
  Directory[]
];
Import[FileNameJoin[{sourceDirectory, "..", "SDPB.m"}]];

(* Keep the spectral data exact until the final numerical conversion. *)
m1 = 1/2;
mgap = 83/50;
mA = 1/5;

J1 = 0;
J2 = 2;

NullIndices = Range[0, 19];
nulllist = {19, -1, -1, -1};
list0 = Table[0, {i, 1, Total[nulllist]+Length[nulllist]}];
functionalDimension = 2 + Length[list0];


Nlist[n_, x_, J_, mA_] := {((-4 mA^2+x)^(3/2) (32 mA^4+2 (-1+J) (8+J) mA^2 x-(-2+J (7+J)) x^2))/(2 x^(5/2)),1/(20 x^(7/2))Sqrt[-4 mA^2+x] (-1280 mA^6+960 mA^4 x+2 (-120+(-1+J) J (7+J) (8+J)) mA^2 x^2-(-20+J (7+J) (-13+J (7+J))) x^3),1/(360 x^(9/2) Sqrt[-4 mA^2+x])(92160 mA^8-92160 mA^6 x+34560 mA^4 x^2+2 (-2880+(-2+J) (-1+J) J (7+J) (8+J) (9+J)) mA^2 x^3-(-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) x^4),1/(10080 x^(11/2) (-4 mA^2+x)^(3/2))(-10321920 mA^10+12902400 mA^8 x-6451200 mA^6 x^2+1612800 mA^4 x^3+2 (-100800+(-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J)) mA^2 x^4+(10080-(-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J)))) x^5),1/(180 x^(9/2) (-4 mA^2+x)^(3/2))J (7+J) (11520 mA^8-11520 mA^6 x-4 (-936+J (7+J) (-26+J (7+J))) mA^4 x^2+2 (-216+J (7+J) (-26+J (7+J))) mA^2 x^3-9 (-23+J (7+J)) x^4),1/(403200 x^(13/2) (-4 mA^2+x)^(5/2))(1651507200 mA^12-2477260800 mA^10 x+1548288000 mA^8 x^2-516096000 mA^6 x^3+96768000 mA^4 x^4+2 (-4838400+(-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J)) mA^2 x^5-(-403200+(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-34+J (5+J)) (-20+J (9+J))) x^6),1/(2520 x^(11/2) (-4 mA^2+x)^(5/2))J (7+J) (-645120 mA^10+806400 mA^8 x-403200 mA^6 x^2-2 (-54720+J (7+J) (924+J (7+J) (-56+J (7+J)))) mA^4 x^3+(-16920+J (7+J) (924+J (7+J) (-56+J (7+J)))) mA^2 x^4-7 (540+J (7+J) (-53+J (7+J))) x^5),1/(360 x^(9/2) (-4 mA^2+x)^(5/2))J (7+J) (2304 (-1+J) (8+J) mA^8+32 (-1+J) (8+J) (-90+J (7+J)) mA^6 x-24 (-1+J) (8+J) (-54+J (7+J)) mA^4 x^2+6 (-1+J) (8+J) (-42+J (7+J)) mA^2 x^3-(540+J (7+J) (-53+J (7+J))) x^4),(-356725555200 mA^14+624269721600 mA^12 x-468202291200 mA^10 x^2+195084288000 mA^8 x^3-48771072000 mA^6 x^4+7315660800 mA^4 x^5+2 (-304819200+(-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J)) mA^2 x^6-(-21772800+(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-62+J (7+J)) (-48+J (7+J)) (-15+J (7+J))) x^7)/(21772800 x^(15/2) (-4 mA^2+x)^(7/2)),1/(100800 x^(13/2) (-4 mA^2+x)^(7/2))J (7+J) (103219200 mA^12-154828800 mA^10 x+96768000 mA^8 x^2-32256000 mA^6 x^3-2 (-2833920+J (7+J) (-44976+J (7+J) (3388+J (7+J) (-100+J (7+J))))) mA^4 x^4+(1+J) (6+J) (-69120+J (7+J) (4024+J (7+J) (-106+J (7+J)))) mA^2 x^5-10 (-31032+J (7+J) (3024+(-7+J) J (7+J) (14+J))) x^6),1/(10080 x^(11/2) (-4 mA^2+x)^(7/2))J (7+J) (-258048 (-1+J) (8+J) mA^10+322560 (-1+J) (8+J) mA^8 x+32 (-1+J) (8+J) (-4500+J (7+J) (-48+J (7+J))) mA^6 x^2-24 (-1+J) (8+J) (-1140+J (7+J) (-48+J (7+J))) mA^4 x^3+6 (-1+J) (8+J) (-300+J (7+J) (-48+J (7+J))) mA^2 x^4-(-31032+J (7+J) (3024+(-7+J) J (7+J) (14+J))) x^5),(99883155456000 mA^16-199766310912000 mA^14 x+174795522048000 mA^12 x^2-87397761024000 mA^10 x^3+27311800320000 mA^8 x^4-5462360064000 mA^6 x^5+682795008000 mA^4 x^6+2 (-24385536000+(-6+J) (-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J) (13+J)) mA^2 x^7-(-1524096000+(-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (-50+J (7+J)) (960+J (7+J) (-83+J (7+J)))) x^8)/(1524096000 x^(17/2) (-4 mA^2+x)^(9/2)),(J (7+J) (-44590694400 mA^14+78033715200 mA^12 x-58525286400 mA^10 x^2+24385536000 mA^8 x^3-6096384000 mA^6 x^4-4 (-240019200+J (7+J) (2888640+J (7+J) (-248256+J (7+J) (9388+J (7+J) (-160+J (7+J)))))) mA^4 x^5+2 (-49507200+J (7+J) (2888640+J (7+J) (-248256+J (7+J) (9388+J (7+J) (-160+J (7+J)))))) mA^2 x^6-27 (1578240+J (7+J) (-209056+J (7+J) (8988+J (7+J) (-160+J (7+J))))) x^7))/(10886400 x^(15/2) (-4 mA^2+x)^(9/2)),-1/(2520 x^(11/2) (-4 mA^2+x)^(9/2))(-2+J) J (7+J) (9+J) (3584 (-1+J) (8+J) mA^10+32 (-10+J) (-1+J) (8+J) (17+J) mA^8 x-32 (-1+J) (8+J) (-100+J (7+J)) mA^6 x^2+4 (-1+J) (8+J) (-230+3 J (7+J)) mA^4 x^3-2 (-1+J) (8+J) (-65+J (7+J)) mA^2 x^4+7 (-53+J (7+J)) x^5),(-35158870720512000 mA^18+79107459121152000 mA^16 x-79107459121152000 mA^14 x^2+46146017820672000 mA^12 x^3-17304756682752000 mA^10 x^4+4326189170688000 mA^8 x^5-721031528448000 mA^6 x^6+77253378048000 mA^4 x^7+2 (-2414168064000+(-7+J) (-6+J) (-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J) (13+J) (14+J)) mA^2 x^8-(-134120448000+(-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (5001600+J (7+J) (-600160+J (7+J) (19516+J (7+J) (-240+J (7+J)))))) x^9)/(134120448000 x^(19/2) (-4 mA^2+x)^(11/2)),(J (7+J) (12485394432000 mA^16-24970788864000 mA^14 x+21849440256000 mA^12 x^2-10924720128000 mA^10 x^3+3413975040000 mA^8 x^4-682795008000 mA^6 x^5-4 (-20447769600+J (7+J) (-236718720+J (7+J) (22252608+J (7+J) (-980520+J (7+J) (21868+J (7+J) (-238+J (7+J))))))) mA^4 x^6+2 (-2158617600+J (7+J) (-236718720+J (7+J) (22252608+J (7+J) (-980520+J (7+J) (21868+J (7+J) (-238+J (7+J))))))) mA^2 x^7-35 (-131466240+J (7+J) (17720496+J (7+J) (-915804+J (7+J) (21808+J (7+J) (-241+J (7+J)))))) x^8))/(762048000 x^(17/2) (-4 mA^2+x)^(11/2)),1/(50400 x^(13/2) (-4 mA^2+x)^(11/2))(-2+J) J (7+J) (9+J) (286720 (-1+J) (8+J) mA^12-430080 (-1+J) (8+J) mA^10 x-16 (-1+J) (8+J) (-15480+J (7+J) (-74+J (7+J))) mA^8 x^2+16 (-1+J) (8+J) (-4280+J (7+J) (-74+J (7+J))) mA^6 x^3-6 (-1+J) (8+J) (-1480+J (7+J) (-74+J (7+J))) mA^4 x^4+(-1+J) (8+J) (-360+J (7+J) (-74+J (7+J))) mA^2 x^5-10 (2564+J (7+J) (-108+J (7+J))) x^6),(15188632151261184000 mA^20-37971580378152960000 mA^18 x+42718027925422080000 mA^16 x^2-28478685283614720000 mA^14 x^3+12459424811581440000 mA^12 x^4-3737827443474432000 mA^10 x^5+778714050723840000 mA^8 x^6-111244864389120000 mA^6 x^7+10429206036480000 mA^4 x^8+2 (-289700167680000+(-8+J) (-7+J) (-6+J) (-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J) (13+J) (14+J) (15+J)) mA^2 x^9-(-14485008384000+(-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (5277600+J (7+J) (-651672+J (7+J) (20980+J (7+J) (-250+J (7+J)))))) x^10)/(14485008384000 x^(21/2) (-4 mA^2+x)^(13/2)),(J (7+J) (-2197429420032000 mA^18+4944216195072000 mA^16 x-4944216195072000 mA^14 x^2+2884126113792000 mA^12 x^3-1081547292672000 mA^10 x^4+270386823168000 mA^8 x^5-45064470528000 mA^6 x^6-2 (-2501346355200+J (7+J) (24088008960+J (7+J) (-2417474304+J (7+J) (118343568+J (7+J) (-3123584+J (7+J) (45192+J (7+J) (-336+J (7+J)))))))) mA^4 x^7+(-388949299200+J (7+J) (24088008960+J (7+J) (-2417474304+J (7+J) (118343568+J (7+J) (-3123584+J (7+J) (45192+J (7+J) (-336+J (7+J)))))))) mA^2 x^8-22 (11405836800+J (7+J) (-1826254080+J (7+J) (107801568+J (7+J) (-3100260+J (7+J) (46228+J (7+J) (-343+J (7+J))))))) x^9))/(33530112000 x^(19/2) (-4 mA^2+x)^(13/2)),((-2+J) J (7+J) (9+J) (-123863040 (-1+J) (8+J) mA^14+216760320 (-1+J) (8+J) mA^12 x-162570240 (-1+J) (8+J) mA^10 x^2-32 (-1+J) (8+J) (-2196000+J (7+J) (5760+J (7+J) (-134+J (7+J)))) mA^8 x^3+32 (-1+J) (8+J) (-608400+J (7+J) (5760+J (7+J) (-134+J (7+J)))) mA^6 x^4-12 (-1+J) (8+J) (-290880+J (7+J) (5760+J (7+J) (-134+J (7+J)))) mA^4 x^5+2 (-1+J) (8+J) (-185040+J (7+J) (5760+J (7+J) (-134+J (7+J)))) mA^2 x^6-27 (-216000+J (7+J) (10752+J (7+J) (-182+J (7+J)))) x^7))/(5443200 x^(15/2) (-4 mA^2+x)^(13/2))}[[n+1]];


(* ---------------------------------------------------------------------- *)
(* PMP assembly                                                           *)
(* ---------------------------------------------------------------------- *)

(* Dense low spins locate the usual extremal states.  The sparse probes
   monitor the transition to the analytic asymptotic blocks.  A final
   solution must still be checked for omitted finite-spin violations. *)
Jmax = 200;
JList = {250, 300, 400, 500, 700, 1000, 1500, 2500, 5000};
continuumJList = DeleteDuplicates @ Join[
  Range[0, Jmax, 2],
  JList
];

(* Positive blockwise rescaling for the q^4 = O(J^8) growth, where
   q = J (J + 7) is the ten-dimensional spin Casimir. *)
spinScale[J_Integer] := 1/(1 + J (J + 7))^4;

(* Analytic boundaries of the large-spin spectrum.  The fixed-mass
   vector was divided by its strictly positive mass-dependent factor. *)
fixedMassLargeJVector = {0, 0, 0, 0, 0, 0, -1};

impactVector[r_] := {
  2,
  0,
  -r/2,
  -r^2/20,
  -r^3/360,
  0,
  -r^4/10080
};

LaunchKernels[];

PMP2SDP[datfile_, prec_: 600] := Module[
    {
      continuumBlocks, specialBlocks, asymptoticBlocks,
      pols, norm, obj, expectedBlocks
    },

    Poly[j_, x_, y_] := Module[{g0, lambda22, polys},
        g0 = 2 x^5 (-4 mA^2 + x)^5;

        (* Only the isolated second state contributes to lambda22. *)
        lambda22 = 0;

        polys = Join[
            {g0, lambda22},
            Nlist[#, x, j, mA] & /@ NullIndices
        ];

        PositiveMatrixWithPrefactor[
            DampedRational[1, {}, 1/E, y],
            {{spinScale[j] polys}}
        ]
    ];

    (* The first isolated state. *)
    Poly1st[j_, x_, y_] := Module[{g0, lambda22, polys},
        g0 = 2 x^5 (-4 mA^2 + x)^5;
        lambda22 = 0;

        polys = Join[
            {g0, lambda22},
            Nlist[#, x, j, mA] & /@ NullIndices
        ];

        PositiveMatrixWithPrefactor[
            DampedRational[1, {}, 1/E, y],
            {{polys}}
        ]
    ];


    (* The isolated spin-2 state. *)
    Poly2nd[j_, x_, y_] := Module[{g0, lambdaH2, polys},
        g0 = 2 x^5 (-4 mA^2 + x)^5;
        lambdaH2 = 1;

        polys = Join[
            {g0, lambdaH2},
            Nlist[#, x, j, mA] & /@ NullIndices
        ];

        PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, y],
            {{polys}}
        ]
    ];

    (* Fixed mass with J -> Infinity. *)
    FixedMassLargeJBlock[y_] := PositiveMatrixWithPrefactor[
      DampedRational[1, {}, 1/E, y],
      {{fixedMassLargeJVector}}
    ];

    (* Joint limit with r = J (J + 7)/s fixed. *)
    ImpactBlock[y_] := PositiveMatrixWithPrefactor[
      DampedRational[1, {}, 1/E, y],
      {{impactVector[y]}}
    ];

    specialBlocks = {
      Poly1st[J1, m1, x],
      Poly2nd[J2, 1, x]
    };

    (* One polynomial block per spin covers the complete mass interval
       s = mgap + x, x >= 0; no mass sampling is required. *)
    continuumBlocks = ParallelMap[
      Poly[#, mgap + x, x] &,
      continuumJList
    ];

    asymptoticBlocks = {
      FixedMassLargeJBlock[x],
      ImpactBlock[x]
    };

    pols = N[
      Join[specialBlocks, continuumBlocks, asymptoticBlocks],
      prec
    ];

    expectedBlocks = 4 + Length[continuumJList];
    If[Length[pols] =!= expectedBlocks,
      Print[
        "Unexpected block count: built ", Length[pols],
        ", expected ", expectedBlocks, "."
      ];
      Abort[]
    ];

    norm = -N[Flatten[{{0, 1}, list0}], prec];
    obj = -N[Flatten[{{1, 0}, list0}], prec];

    If[
        Length[norm] =!= functionalDimension ||
        Length[obj] =!= functionalDimension,
        Print["Objective or normalization dimension mismatch."];
        Abort[]
    ];

    Print["functional dimension = ", functionalDimension];
    Print["number of PMP blocks = ", Length[pols] ];
    Print["finite J values = ", continuumJList];
    Print["Included fixed-mass and correlated large-J blocks."];

    WritePmpJson[
        datfile,
        SDP[obj, norm, pols],
        prec,
        getAnalyticSampleData
    ]
];

outputFile = FileNameJoin[{sourceDirectory, "n_pmp.json"}];
If[!TrueQ[$Test16SkipExport],
  PMP2SDP[outputFile, 600]
];
