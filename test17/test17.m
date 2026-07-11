Import["../SDPB.m"]

ClearAll[NumericalPositiveMatrixWithPrefactor];

toJsonNestedNumberArray[expr_, prec_] :=
  expr /. n_?NumericQ :> toJsonNumber[n, prec];

toJsonObject[
  NumericalPositiveMatrixWithPrefactor[pmp_?AssociationQ],
  prec_,
  getSampleDataFn_:Function[<||>]
] := Module[
  {sampleData, samplePoints, sampleScalings, functionValues, basisValues,
   prefactor, reducedPrefactor},
  sampleData = getSampleDataFn[NumericalPositiveMatrixWithPrefactor[pmp], prec];
  prefactor = Lookup[pmp, "prefactor", 1];
  reducedPrefactor = Lookup[pmp, "reducedPrefactor", prefactor];
  samplePoints = Lookup[sampleData, "samplePoints",
    Lookup[pmp, "samplePoints", Missing[]]];
  sampleScalings = Lookup[sampleData, "sampleScalings",
    Lookup[pmp, "sampleScalings", Missing[]]];
  functionValues = Lookup[pmp, "polynomials", Missing[]];
  basisValues = Lookup[sampleData, "basisValues",
    Lookup[pmp, "basisValues", Missing[]]];
  DeleteMissing @ <|
    "prefactor" -> toJsonDampedRational[prefactor, prec],
    "reducedPrefactor" -> toJsonDampedRational[reducedPrefactor, prec],
    "samplePoints" -> toJsonNumberArray[samplePoints, prec],
    "sampleScalings" -> toJsonNumberArray[sampleScalings, prec],
    "polynomials" -> If[MissingQ[functionValues], Missing[],
      toJsonNestedNumberArray[functionValues, prec]],
    "basisValues" -> If[MissingQ[basisValues], Missing[],
      toJsonNestedNumberArray[basisValues, prec]]
  |>
];

WritePmpJsonNumerical[
  file_,
  SDP[objective_, normalization_, positiveMatricesWithPrefactors_],
  prec_,
  getSampleDataFn_:Function[<||>]
] := exportJson[
  file,
  <|
    "objective" -> toJsonNumberArray[objective, prec],
    "normalization" -> toJsonNumberArray[normalization, prec],
    "PositiveMatrixWithPrefactorArray" ->
      Table[toJsonObject[pmp, prec, getSampleDataFn],
        {pmp, positiveMatricesWithPrefactors}]
  |>
];
m1 = N[2/5, 1000];
J1 = 0;
J2 = 2;
mgap = N[166/100, 1000];

nulllist = {51, -1, -1, -1};
list0 = Table[0, {i, 1, Total[nulllist]+Length[nulllist]}];

Nlist[n_, z_, J_] := {(2-J (7+J))/(2 z^2),(20-J (7+J) (-13+J (7+J)))/(20 z^3),-(((-12+J (7+J)) (30+J (7+J) (-23+J (7+J))))/(360 z^4)),1/z^5-((-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J))))/(10080 z^5),-((J (7+J) (-23+J (7+J)))/(20 z^5)),1/z^6-((-4+J) (-2+J) J (7+J) (9+J) (11+J) (-34+J (5+J)) (-20+J (9+J)))/(403200 z^6),-((J (7+J) (540+J (7+J) (-53+J (7+J))))/(360 z^6)),1/z^7-((-4+J) (-2+J) J (7+J) (9+J) (11+J) (-62+J (7+J)) (-48+J (7+J)) (-15+J (7+J)))/(21772800 z^7),-((J (7+J) (-31032+J (7+J) (3024+(-7+J) J (7+J) (14+J))))/(10080 z^7)),1/z^8-((-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (-50+J (7+J)) (960+J (7+J) (-83+J (7+J))))/(1524096000 z^8),-((J (7+J) (1578240+J (7+J) (-209056+J (7+J) (8988+J (7+J) (-160+J (7+J))))))/(403200 z^8)),-(((-2+J) J (7+J) (9+J) (-53+J (7+J)))/(360 z^8)),1/z^9-1/(134120448000 z^9) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (5001600+J (7+J) (-600160+J (7+J) (19516+J (7+J) (-240+J (7+J))))),-((J (7+J) (-131466240+J (7+J) (17720496+J (7+J) (-915804+J (7+J) (21808+J (7+J) (-241+J (7+J)))))))/(21772800 z^9)),-(((-2+J) J (7+J) (9+J) (2564+J (7+J) (-108+J (7+J))))/(5040 z^9)),1/z^10-1/(14485008384000 z^10) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (5277600+J (7+J) (-651672+J (7+J) (20980+J (7+J) (-250+J (7+J))))),-(1/(1524096000 z^10))J (7+J) (11405836800+J (7+J) (-1826254080+J (7+J) (107801568+J (7+J) (-3100260+J (7+J) (46228+J (7+J) (-343+J (7+J))))))),-(((-2+J) J (7+J) (9+J) (-216000+J (7+J) (10752+J (7+J) (-182+J (7+J)))))/(201600 z^10)),1/z^11-1/(1883051089920000 z^11) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (-797572800+J (7+J) (106776000+J (7+J) (-3946556+J (7+J) (60108+J (7+J) (-405+J (7+J)))))),-(1/(134120448000 z^11))J (7+J) (-1379752704000+J (7+J) (225934848000+J (7+J) (-14770082880+J (7+J) (486509168+J (7+J) (-8812960+J (7+J) (88928+J (7+J) (-468+J (7+J)))))))),-(((-2+J) J (7+J) (9+J) (21948480+J (7+J) (-1323000+J (7+J) (28702+J (7+J) (-277+J (7+J))))))/(10886400 z^11)),-(((-2+J) J (7+J) (9+J) (35064000+J (7+J) (-1763640+J (7+J) (31942+J (7+J) (-277+J (7+J))))))/(21772800 z^11)),1/z^12-1/(289989867847680000 z^12) (-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (-833676480+J (7+J) (114378912+J (7+J) (-4230052+J (7+J) (63448+J (7+J) (-417+J (7+J)))))),-(1/(14485008384000 z^12))J (7+J) (180894269952000+J (7+J) (-33122195274240+J (7+J) (2349433112832+J (7+J) (-85849259040+J (7+J) (1788748560+J (7+J) (-22054832+J (7+J) (158952+J (7+J) (-618+J (7+J))))))))),-(1/(762048000 z^12))(-2+J) J (7+J) (9+J) (-2563473600+J (7+J) (175893840+J (7+J) (-4726576+J (7+J) (61658+J (7+J) (-395+J (7+J)))))),-(1/(1524096000 z^12))(-2+J) J (7+J) (9+J) (-5573865600+J (7+J) (318324240+J (7+J) (-6824476+J (7+J) (71108+J (7+J) (-395+J (7+J)))))),1/z^13-1/(52198176212582400000 z^13) (-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (172008748800+J (7+J) (-25008497280+J (7+J) (1017235728+J (7+J) (-17785440+J (7+J) (152128+J (7+J) (-628+J (7+J))))))),-(1/(1883051089920000 z^13))J (7+J) (-30370283513856000+J (7+J) (5685348310272000+J (7+J) (-431278994177280+J (7+J) (17117281051776+J (7+J) (-396920087856+J (7+J) (5654064144+J (7+J) (-50058872+J (7+J) (268176+J (7+J) (-795+J (7+J)))))))))),-(1/(67060224000 z^13))(-2+J) J (7+J) (9+J) (357678604800+J (7+J) (-27684313920+J (7+J) (855078608+J (7+J) (-13637120+J (7+J) (118668+J (7+J) (-538+J (7+J))))))),-(1/(134120448000 z^13))(-2+J) J (7+J) (9+J) (942637132800+J (7+J) (-59261829120+J (7+J) (1465072608+J (7+J) (-18734520+J (7+J) (134068+J (7+J) (-538+J (7+J))))))),1/z^14-((-12+J) (-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (19+J) (178803072000+J (7+J) (-26522288640+J (7+J) (1083860832+J (7+J) (-18822392+J (7+J) (158628+J (7+J) (-642+J (7+J))))))))/(10857220652217139200000 z^14),-(1/(289989867847680000 z^14))J (7+J) (5582653171070976000+J (7+J) (-1130338166370048000+J (7+J) (90722029433978880+J (7+J) (-3853417628059776+J (7+J) (97298039791872+J (7+J) (-1547871265360+J (7+J) (15899813776+J (7+J) (-105145568+J (7+J) (431816+J (7+J) (-1001+J (7+J))))))))))),-(1/(289989867847680000 z^14))J (7+J) (18904909351071744000+J (7+J) (-3159997781508403200+J (7+J) (212491113563151360+J (7+J) (-7648220048749056+J (7+J) (164888697966912+J (7+J) (-2260267109520+J (7+J) (20284193776+J (7+J) (-119680088+J (7+J) (451836+J (7+J) (-1001+J (7+J))))))))))),-(1/(7242504192000 z^14))(-2+J) J (7+J) (9+J) (-57858694502400+J (7+J) (4940899119360+J (7+J) (-172748862720+J (7+J) (3196589328+J (7+J) (-34081280+J (7+J) (211008+J (7+J) (-708+J (7+J)))))))),-(1/(14485008384000 z^14))(-2+J) J (7+J) (9+J) (-179566979328000+J (7+J) (12489984864000+J (7+J) (-350374821120+J (7+J) (5203121328+J (7+J) (-45129680+J (7+J) (234768+J (7+J) (-708+J (7+J)))))))),1/z^15-((-12+J) (-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (19+J) (-48228426086400+J (7+J) (7471309098240+J (7+J) (-327436493664+J (7+J) (6325658024+J (7+J) (-62930772+J (7+J) (336322+J (7+J) (-917+J (7+J)))))))))/(2584018515227679129600000 z^15),-(1/(52198176212582400000 z^15))J (7+J) (-1247747313202790400000+J (7+J) (257807412952510464000+J (7+J) (-21713822416662681600+J (7+J) (976289459242567680+J (7+J) (-26434504759413888+J (7+J) (459223576873344+J (7+J) (-5286173524128+J (7+J) (40717222480+J (7+J) (-207319640+J (7+J) (668976+J (7+J) (-1238+J (7+J)))))))))))),-(1/(52198176212582400000 z^15))J (7+J) (-4427589781962547200000+J (7+J) (778134090027525120000+J (7+J) (-55587311278507929600+J (7+J) (2133844925721868800+J (7+J) (-49603448066289408+J (7+J) (744185486893824+J (7+J) (-7462993855968+J (7+J) (50767274800+J (7+J) (-232960640+J (7+J) (696696+J (7+J) (-1238+J (7+J)))))))))))),-(1/(941525544960000 z^15))(-2+J) J (7+J) (9+J) (10916797731840000+J (7+J) (-1018511499110400+J (7+J) (39142949166720+J (7+J) (-814016017152+J (7+J) (10077356168+J (7+J) (-76691252+J (7+J) (353250+J (7+J) (-907+J (7+J))))))))),-(1/(1883051089920000 z^15))(-2+J) J (7+J) (9+J) (38408970023424000+J (7+J) (-2914947731980800+J (7+J) (91204340699520+J (7+J) (-1546371339552+J (7+J) (15657273368+J (7+J) (-98663852+J (7+J) (388350+J (7+J) (-907+J (7+J))))))))),1/z^16-((-14+J) (-12+J) (-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (19+J) (21+J) (-49949104896000+J (7+J) (7864546348800+J (7+J) (-346775577120+J (7+J) (6684386328+J (7+J) (-65927780+J (7+J) (347718+J (7+J) (-933+J (7+J)))))))))/(697684999111473364992000000 z^16),-((J (7+J) (305478770084167680000000+J (7+J) (-66880855243136286720000+J (7+J) (5873152091579080704000+J (7+J) (-277205314388051865600+J (7+J) (7958521963421475840+J (7+J) (-148608950149018368+J (7+J) (1873341464128704+J (7+J) (-16222546312128+J (7+J) (96556311280+J (7+J) (-387804560+J (7+J) (1003236+J (7+J) (-1508+J (7+J))))))))))))))/(10857220652217139200000 z^16)),-((J (7+J) (1188386342061144145920000+J (7+J) (-217971526056244789248000+J (7+J) (16315470556851061555200+J (7+J) (-662248807410304512000+J (7+J) (16420500781059271680+J (7+J) (-265833327608356608+J (7+J) (2921409120924864+J (7+J) (-22249607260608+J (7+J) (118056830320+J (7+J) (-431047760+J (7+J) (1040676+J (7+J) (-1508+J (7+J))))))))))))))/(10857220652217139200000 z^16)),-(1/(144994933923840000 z^16))(-2+J) J (7+J) (9+J) (-2358577574825472000+J (7+J) (237610799759116800+J (7+J) (-9941557443736320+J (7+J) (227429331438240+J (7+J) (-3164274187792+J (7+J) (28016444448+J (7+J) (-159178952+J (7+J) (563810+J (7+J) (-1137+J (7+J)))))))))),-(1/(289989867847680000 z^16))(-2+J) J (7+J) (9+J) (-9288820773872640000+J (7+J) (765214719981696000+J (7+J) (-26297097573945600+J (7+J) (496256191980480+J (7+J) (-5707230603792+J (7+J) (41913607728+J (7+J) (-200019752+J (7+J) (613860+J (7+J) (-1137+J (7+J)))))))))),1/z^17-((-14+J) (-12+J) (-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (19+J) (21+J) (17046656188416000+J (7+J) (-2775655413504000+J (7+J) (129135434676480+J (7+J) (-2693784807936+J (7+J) (29839696528+J (7+J) (-187759616+J (7+J) (673224+J (7+J) (-1280+J (7+J))))))))))/(212096239729887902957568000000 z^17),-((J (7+J) (-87853552280843169300480000+J (7+J) (19584641136608558678016000+J (7+J) (-1783834766769185877196800+J (7+J) (87808102998808260096000+J (7+J) (-2650174153472579208192+J (7+J) (52576896713946183936+J (7+J) (-714113670487896000+J (7+J) (6790569461809664+J (7+J) (-45576671260848+J (7+J) (214685735264+J (7+J) (-693742452+J (7+J) (1463280+J (7+J) (-1813+J (7+J)))))))))))))))/(2584018515227679129600000 z^17)),-((J (7+J) (-356782416111642673152000000+J (7+J) (68177518035168076333056000+J (7+J) (-5344027626418386709708800+J (7+J) (227778687728448460369920+J (7+J) (-5972738835680399729664+J (7+J) (103209702498806706432+J (7+J) (-1225240506943814592+J (7+J) (10240675349197824+J (7+J) (-60956946250928+J (7+J) (258094604768+J (7+J) (-763939124+J (7+J) (1512784+J (7+J) (-1813+J (7+J)))))))))))))))/(2584018515227679129600000 z^17)),-(1/(26099088106291200000 z^17))(-2+J) J (7+J) (9+J) (582699698380800000000+J (7+J) (-62900634396722688000+J (7+J) (2822389999565068800+J (7+J) (-69946031480532480+J (7+J) (1070137386647136+J (7+J) (-10661654391696+J (7+J) (70646061304+J (7+J) (-309730172+J (7+J) (865536+J (7+J) (-1400+J (7+J))))))))))),-(1/(52198176212582400000 z^17))(-2+J) J (7+J) (9+J) (2521949401939968000000+J (7+J) (-223540198278644736000+J (7+J) (8341871152470912000+J (7+J) (-172808387428884480+J (7+J) (2211002533791936+J (7+J) (-18416939498496+J (7+J) (102380333104+J (7+J) (-381594272+J (7+J) (934836+J (7+J) (-1400+J (7+J))))))))))),-(1/(627683696640000 z^17))(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-1098365329152000+J (7+J) (57080097100800+J (7+J) (-1240016055840+J (7+J) (14651007048+J (7+J) (-102462260+J (7+J) (427338+J (7+J) (-993+J (7+J)))))))),1/z^18-((-16+J) (-14+J) (-12+J) (-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (19+J) (21+J) (23+J) (17607094156800000+J (7+J) (-2905063330752000+J (7+J) (136020135744000+J (7+J) (-2838087378720+J (7+J) (31287521168+J (7+J) (-195162000+J (7+J) (691768+J (7+J) (-1298+J (7+J))))))))))/(72112721508161887005573120000000 z^18)}[[n+1]];


MNlist[n_,z_,J_] := {
	{0,0,0},
	{0,Nlist[n,z,J],0},
	{0,0,0}
};


polyify[expr_] := Expand @ Cancel @ Together[expr];

NPolyInf[n_,J_,x_] := {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-(1/72112721508161887005573120000000)}[[n+1]];


LaunchKernels[];

ClearAll[
  FlattenSpinTiers, UToZ, ZToU, RegularScalarVector, ScalarVectorToMatrix,
  AnalyticPoly, AnalyticPolyInf, AnalyticPoly2nd, NumericalScalarBlock,
  ChebyshevLobattoNodes, ToInterval, ChebyshevCoefficients,
  ChebyshevApproxValue, ValidationNodes, NormalizeVector,
  EstimateComponentScales, AnalyzeIntervalForSpin, HarvestIntervalPoints,
  AnalyzeNonPolynomialFamily, BuildCompressedSamplingGrid,
  BlocksFromSamplingGrid, FunctionalEntryNonzeroQ
];

FlattenSpinTiers[jTiers_] := Join @@ jTiers;

UToZ[u_] := mgap + Exp[u] - 1;
ZToU[z_] := Log[1 + z - mgap];

RegularScalarVector[j_, z_] := Module[{pref = z^18},
  Join[
    {polyify[pref*(2/z)], polyify[pref*0]},
    Table[polyify[pref*Nlist[n, z, j]], {n, 0, nulllist[[1]]}]
  ]
];

ScalarVectorToMatrix[values_] :=
  Table[
    Table[
      If[row == 2 && column == 2, values, Table[0, {Length[values]}]],
      {column, 3}
    ],
    {row, 3}
  ];

AnalyticPoly[j_, z_, y_] := PositiveMatrixWithPrefactor[
  DampedRational[1, {}, 1/E, y],
  ScalarVectorToMatrix[RegularScalarVector[j, z]]
];

AnalyticPolyInf[j_, z_, y_] := PositiveMatrixWithPrefactor[
  DampedRational[1, {}, 1/E, y],
  ScalarVectorToMatrix[
    Join[
      {0, 0},
      Table[NPolyInf[n, j, z], {n, 0, nulllist[[1]]}]
    ]
  ]
];

AnalyticPoly2nd[j_, z_, y_] := Module[{pref = z^18},
  PositiveMatrixWithPrefactor[
    DampedRational[1, {}, 1/E, y],
    ScalarVectorToMatrix[
      Join[
        {polyify[pref*(2/z)], polyify[pref]},
        Table[polyify[pref*Nlist[n, z, j]], {n, 0, nulllist[[1]]}]
      ]
    ]
  ]
];

NumericalScalarBlock[j_, u0_, vectorFn_, prec_] := Module[{values},
  values = N[vectorFn[j, UToZ[u0]], prec];
  NumericalPositiveMatrixWithPrefactor[<|
    "prefactor" -> DampedRational[1, {}, 1/E, u0],
    "samplePoints" -> {u0},
    "sampleScalings" -> {Exp[-u0]},
    "polynomials" ->
      Table[
        Table[
          If[row == 2 && column == 2,
            Table[{N[value, prec]}, {value, values}],
            Table[{0}, {Length[values]}]
          ],
          {column, 3}
        ],
        {row, 3}
      ]
  |>]
];

ChebyshevLobattoNodes[p_Integer?Positive] :=
  Reverse[Table[Cos[k Pi/p], {k, 0, p}]];

ToInterval[xi_, {a_, b_}] := (a + b)/2 + (b - a) xi/2;

ChebyshevCoefficients[values_List] := Module[
  {p = Length[values] - 1, coeffs, endpointWeight},
  endpointWeight[j_] := If[j == 0 || j == p, 1/2, 1];
  coeffs = Table[
    (2/p) Total[
      Table[
        endpointWeight[j] values[[j + 1]] Cos[k j Pi/p],
        {j, 0, p}
      ]
    ],
    {k, 0, p}
  ];
  ReplacePart[coeffs, {1 -> coeffs[[1]]/2, Length[coeffs] -> Last[coeffs]/2}]
];

ChebyshevApproxValue[coeffs_List, xi_] :=
  Total[
    Table[coeffs[[k + 1]] ChebyshevT[k, xi], {k, 0, Length[coeffs] - 1}]
  ];

ValidationNodes[p_Integer?Positive] :=
  Table[Cos[(k + 1/2) Pi/(p + 1)], {k, 0, p}];

NormalizeVector[values_, scales_] := N[values/scales];

EstimateComponentScales[spins_, vectorFn_, {uMin_, uMax_}, probeCount_, prec_] :=
 Module[{probeUs, rows, columns},
  probeUs = Subdivide[uMin, uMax, probeCount];
  rows = Flatten[
    Table[
      Abs[N[vectorFn[j, UToZ[uu]], Min[prec, 80]]],
      {j, spins}, {uu, probeUs}
    ],
    1
  ];
  columns = Transpose[rows];
  Max[1, Max[#]] & /@ columns
];

AnalyzeIntervalForSpin[j_, interval_, vectorFn_, scales_, config_] := Module[
  {
    orders = Lookup[config, "orders", {8, 16, 24, 32}],
    tailCount = Lookup[config, "tailCount", 4],
    spectralTol = Lookup[config, "spectralTolerance", 10^-14],
    validationTol = Lookup[config, "validationTolerance", 10^-12],
    prec = Lookup[config, "precision", 200],
    accepted = False, result, p, nodes, us, values, coeffs,
    tailNorm, totalNorm, spectralError, xis, vUs, actual, interp,
    validationErrors, validationError
  },
  Do[
    nodes = ChebyshevLobattoNodes[p];
    us = ToInterval[#, interval] & /@ nodes;
    values = NormalizeVector[N[vectorFn[j, UToZ[#]], Min[prec, 80]], scales] & /@ us;
    coeffs = ChebyshevCoefficients[values];
    tailNorm = Max[
      Total[
        Abs[Take[coeffs, -Min[tailCount, Length[coeffs]]]]
      ]
    ];
    totalNorm = Max[1, Max[Total[Abs[coeffs]]]];
    spectralError = N[tailNorm/totalNorm];

    xis = ValidationNodes[p];
    vUs = ToInterval[#, interval] & /@ xis;
    validationErrors = Table[
      actual = NormalizeVector[N[vectorFn[j, UToZ[vUs[[i]]]], Min[prec, 80]], scales];
      interp = ChebyshevApproxValue[coeffs, xis[[i]]];
      Norm[actual - interp, Infinity]/Max[1, Norm[actual, Infinity]],
      {i, Length[xis]}
    ];
    validationError = Max[validationErrors];
    result = <|
      "spin" -> j,
      "interval" -> interval,
      "order" -> p,
      "nodes" -> us,
      "coefficients" -> coeffs,
      "validationNodes" -> vUs,
      "validationErrors" -> validationErrors,
      "spectralError" -> spectralError,
      "validationError" -> validationError,
      "accepted" -> TrueQ[
        spectralError <= spectralTol && validationError <= validationTol
      ]
    |>;
    If[result["accepted"], accepted = True; Break[]],
    {p, orders}
  ];
  result
];

HarvestIntervalPoints[intervalData_, config_] := Module[
  {
    interval = intervalData["interval"],
    guardFraction = Lookup[config, "guardFraction", 0.03],
    keepAllNodes = Lookup[config, "keepChebyshevNodes", True],
    nodes, errors, maxErrorU, h, points
  },
  h = interval[[2]] - interval[[1]];
  nodes = If[keepAllNodes, intervalData["nodes"], {}];
  errors = intervalData["validationErrors"];
  maxErrorU = If[Length[errors] > 0,
    intervalData["validationNodes"][[First@Ordering[errors, -1]]],
    Mean[interval]
  ];
  points = Join[
    interval,
    nodes,
    {maxErrorU - guardFraction h, maxErrorU, maxErrorU + guardFraction h}
  ];
  Select[points, interval[[1]] <= # <= interval[[2]] &]
];

AnalyzeNonPolynomialFamily[config_] := Module[
  {
    families = Lookup[config, "families", {}],
    spins = Lookup[config, "spins", {}],
    uRange = Lookup[config, "uRange", {0, Log[200]}],
    initialIntervals = Lookup[config, "initialIntervals", 8],
    maxIntervals = Lookup[config, "maxIntervalsPerSpin", 64],
    prec = Lookup[config, "precision", 200],
    probeCount = Lookup[config, "scaleProbeCount", 32],
    tol = Lookup[config, "intervalTolerance", 10^-12],
    family, vectorFn, scales, intervalsBySpin, dataBySpin,
    samplesBySpin, unresolvedBySpin, intervals, queue, current, analyzed,
    halves, retainedData, grid
  },
  If[families === {} || spins === {},
    Return[<|
      "scales" -> {},
      "reducedBasis" -> <|"rank" -> 0, "singularValues" -> {}|>,
      "spinClusters" -> <||>,
      "localApproximants" -> <||>,
      "criticalPoints" -> <||>,
      "tailData" -> <|"uRange" -> uRange, "analyticTailCertified" -> False|>,
      "denseProbeDiagnostics" -> <|"message" -> "No genuinely non-polynomial families were configured."|>,
      "samplesBySpin" -> <||>
    |>]
  ];

  family = First[families];
  vectorFn = family["function"];
  scales = EstimateComponentScales[spins, vectorFn, uRange, probeCount, prec];

  intervalsBySpin = Association[];
  dataBySpin = Association[];
  samplesBySpin = Association[];
  unresolvedBySpin = Association[];

  Do[
    queue = Partition[Subdivide[uRange[[1]], uRange[[2]], initialIntervals], 2, 1];
    retainedData = {};
    While[Length[queue] > 0 && Length[retainedData] + Length[queue] <= maxIntervals,
      current = First[queue];
      queue = Rest[queue];
      analyzed = AnalyzeIntervalForSpin[j, current, vectorFn, scales, config];
      If[analyzed["accepted"] || (current[[2]] - current[[1]] <= tol),
        retainedData = Append[retainedData, analyzed],
        halves = {{current[[1]], Mean[current]}, {Mean[current], current[[2]]}};
        queue = Join[halves, queue]
      ];
    ];
    If[Length[queue] > 0,
      AssociateTo[unresolvedBySpin, j -> queue];
      retainedData = Join[
        retainedData,
        AnalyzeIntervalForSpin[j, #, vectorFn, scales, config] & /@ queue
      ];
    ];
    grid = DeleteDuplicates[
      Sort[Flatten[HarvestIntervalPoints[#, config] & /@ retainedData]],
      Abs[#1 - #2] < Lookup[config, "duplicateTolerance", 10^-20] &
    ];
    AssociateTo[intervalsBySpin, j -> retainedData[[All, "interval"]]];
    AssociateTo[dataBySpin, j -> retainedData];
    AssociateTo[samplesBySpin, j -> SetPrecision[grid, prec]],
    {j, spins}
  ];

  <|
    "scales" -> scales,
    "reducedBasis" -> <|
      "rank" -> Length[scales],
      "note" -> "Componentwise normalization is active; SVD truncation is intentionally disabled unless a genuine non-polynomial family is supplied and audited."
    |>,
    "spinClusters" -> AssociationThread[spins, List /@ spins],
    "localApproximants" -> dataBySpin,
    "criticalPoints" -> samplesBySpin,
    "tailData" -> <|
      "uRange" -> uRange,
      "zRange" -> UToZ /@ uRange,
      "analyticTailCertified" -> False,
      "note" -> "The finite tail is extended by the chosen uRange; no rigorous asymptotic remainder certificate is claimed here."
    |>,
    "denseProbeDiagnostics" -> <|
      "unresolvedIntervalsBySpin" -> unresolvedBySpin,
      "maxValidationError" -> Max[
        Flatten[
          Values[dataBySpin][[All, All, "validationError"]]
        ]
      ],
      "maxSpectralError" -> Max[
        Flatten[
          Values[dataBySpin][[All, All, "spectralError"]]
        ]
      ]
    |>,
    "samplesBySpin" -> samplesBySpin
  |>
];

BuildCompressedSamplingGrid[analysis_, config_] := Module[
  {
    samplesBySpin = analysis["samplesBySpin"],
    spins = Lookup[config, "spins", Keys[analysis["samplesBySpin"]]],
    oldCommonGridCount = Lookup[config, "oldCommonGridCount", 242],
    maxGrowth = Lookup[config, "maxBlockGrowthOverOldCartesian", 1],
    blockCount, oldCartesianCount
  },
  blockCount = Total[Length /@ Values[samplesBySpin]];
  oldCartesianCount = Length[spins] oldCommonGridCount;
  If[oldCartesianCount > 0 && blockCount > maxGrowth oldCartesianCount,
    Print[
      "ERROR: adaptive sampling produced ", blockCount,
      " blocks, exceeding the allowed growth over the old Cartesian count ",
      oldCartesianCount, ". Tighten compression or loosen diagnostics first."
    ];
    Quit[1]
  ];
  <|
    "samplesBySpin" -> samplesBySpin,
    "blockCount" -> blockCount,
    "oldCartesianCount" -> oldCartesianCount,
    "validationError" -> Lookup[analysis["denseProbeDiagnostics"], "maxValidationError", 0],
    "tailError" -> Missing["NotCertified"],
    "compressionDiagnostics" -> <|
      "method" -> "per-spin Chebyshev grid with duplicate removal",
      "coneCompression" -> "not applied; scalar cone compression should be enabled only after an NNLS/LP residual audit"
    |>
  |>
];

BlocksFromSamplingGrid[gridData_, vectorFn_, prec_] := Module[
  {rules = Normal[gridData["samplesBySpin"]]},
  If[rules === {}, Return[{}]];
  Flatten[
    Table[
      NumericalScalarBlock[rule[[1]], uu, vectorFn, prec],
      {rule, rules},
      {uu, rule[[2]]}
    ],
    1
  ]
];

FunctionalEntryNonzeroQ[pmp_, k_] := Module[{entries},
  entries = Which[
    Head[pmp] === NumericalPositiveMatrixWithPrefactor,
      Flatten[pmp[[1]]["polynomials"][[All, All, k, All]]],
    Head[pmp] === PositiveMatrixWithPrefactor && AssociationQ[pmp[[1]]],
      Flatten[pmp[[1]]["polynomials"][[All, All, k]]],
    Head[pmp] === PositiveMatrixWithPrefactor,
      Flatten[pmp[[2]][[All, All, k]]],
    True,
      {}
  ];
  AnyTrue[entries, ! TrueQ[Simplify[# == 0]] &]
];

PMP2SDP[datfile_, prec_:600] := Module[
    {
        jTiers, allSpins, phaseAConfig, analysis, gridData,
        analyticRegularBlocks, sampledRemainderBlocks, pols, norm, obj,
        functionalCount, functionalCovered, missingFunctionals
    },
    jTiers = {
      Range[0, 1000, 2],
      Range[1500, 5000, 100],
      Range[6000, 20000, 500],
      Range[20000, 50000, 2000]
    };
    allSpins = FlattenSpinTiers[jTiers];

    phaseAConfig = <|
      "precision" -> prec,
      "spins" -> allSpins,
      (* Add genuine non-polynomial remainders here, never the polynomial-in-z regular family. *)
      "families" -> {},
      "uRange" -> {0, Log[200]},
      "orders" -> {8, 16, 24, 32},
      "tailCount" -> 4,
      "spectralTolerance" -> 10^-14,
      "validationTolerance" -> 10^-12,
      "initialIntervals" -> 8,
      "maxIntervalsPerSpin" -> 64,
      "guardFraction" -> 0.03,
      "oldCommonGridCount" -> 242,
      "maxBlockGrowthOverOldCartesian" -> 1
    |>;

    Print["J samples: ", Total[Length /@ jTiers]];
    Print[
      "Regular family is polynomial in z after multiplying by z^18; ",
      "using analytic PMP blocks instead of the old compactified sampling grid."
    ];
    Print["Phase-A sampler is configured for genuinely non-polynomial remainders."];

    analysis = AnalyzeNonPolynomialFamily[phaseAConfig];
    gridData = BuildCompressedSamplingGrid[analysis, phaseAConfig];
    sampledRemainderBlocks = BlocksFromSamplingGrid[
      gridData,
      RegularScalarVector,
      prec
    ];

    analyticRegularBlocks = N[
      ParallelTable[AnalyticPoly[j, mgap + x, x], {j, allSpins}],
      prec
    ];

    pols = Join[
      {N[AnalyticPoly2nd[J2, 1, x], prec]},
      analyticRegularBlocks,
      sampledRemainderBlocks,
      {
        N[AnalyticPoly[J1, m1, x], prec],
        N[AnalyticPolyInf[0, mgap + x, x], prec]
      }
    ];

    Print["Adaptive numerical remainder blocks: ", Length[sampledRemainderBlocks]];
    Print["Built ", Length[pols], " PMP blocks total."];

    norm = -1 * N[Flatten[{{0, 1}, list0}], prec];
    obj = -1 * N[Flatten[{{1, 0}, list0}], prec];

    functionalCount = Length[norm];
    functionalCovered = Table[
      AnyTrue[pols, FunctionalEntryNonzeroQ[#, k] &],
      {k, functionalCount}
    ];
    missingFunctionals = Flatten @ Position[functionalCovered, False];

    If[missingFunctionals =!= {},
      Print[
        "ERROR: sampled PMP has identically zero functional indices: ",
        missingFunctionals,
        ". Expand the x/J sampling grid before running SDPB."
      ];
      Quit[1]
    ];

    Print["Validated all ", functionalCount,
      " functional directions: none are identically zero."];
    Print["size of norm = ", Length[norm]];
    Print["size of obj = ", Length[obj]];
    Print["Writing ", datfile, "..."];
    WritePmpJsonNumerical[datfile, SDP[obj, norm, pols], prec,
      getAnalyticSampleData];
    Print["Wrote ", datfile, "."]
];

PMP2SDP["n_pmp.json", 1000];
