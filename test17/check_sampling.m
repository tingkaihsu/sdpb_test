(*
  MEMORY-BOUNDED ADAPTIVE X REFINEMENT FOR test17.m

  This script owns a hardcoded copy of test17.m's xSamples, using the
  same xSamples = SetPrecision[...] assignment form.
  It never reads sample points from a file, command line, or environment.

  MANUAL TWO-COPY-PASTE WORKFLOW AFTER EACH REFINEMENT:
    1. Run SDPB with test17.m.
    2. Run this script against SDPB's y.txt.
    3. Copy the emitted xSamples = SetPrecision[{...}, prec]; block
       unchanged into BOTH:
         - test17.m's xSamples assignment, and
         - the xSamples assignment below in this script.
       The two copies must remain synchronized by hand.

  Every interval midpoint is tested against every J sampled by test17.m.
  For a negative interval, score = Abs[F(mid)]*(xb-xa).  This score is a
  greedy criticality proxy, not a provably optimal selection of the most
  critical points.  Highest-score intervals are refined first subject to
  maxTotalPoints, with no partial-interval refinement.

  Original points are always merged into the result.  Never discard them:
  doing so can make negative regions oscillate between refinement rounds.

  Usage:
    wolframscript -file refine_sampling_test17.m \
      <yFile> [nPts] [minWidth] <maxTotalPoints>

  Arguments:
    yFile           required path to SDPB y.txt
    nPts            subdivisions per chosen interval (default 10);
                    each interval contributes nPts-1 interior points
    minWidth        unconditional minimum refinable width (default 1e-6)
    maxTotalPoints  required hard cap on the merged, deduplicated grid

  Output file (fixed): refined_xSamples_test17.txt
  Output points are exact Mathematica rationals (for example 1/200),
  not decimal approximations.

  Exit codes:
    0  true convergence: zero refinable negative intervals
    1  at least one interval refined within budget
    2  argument, input, parsing, or evaluation error
    3  refinable negative intervals remain, but budget prevented refinement
*)

prec = 1000;

xSamples = SetPrecision[{0, 1/200, 1/100, 3/200, 31/2000, 2/125, 33/2000, 17/1000, 7/400, 9/500, 37/2000, 19/1000, 39/2000, 1/50, 1/40, 3/100, 61/2000, 31/1000, 63/2000, 4/125, 13/400, 33/1000, 67/2000, 17/500, 69/2000, 7/200, 1/25, 9/200, 1/20, 11/200, 3/50, 121/2000, 61/1000, 123/2000, 31/500, 1/16, 63/1000, 127/2000, 8/125, 129/2000, 13/200, 7/100, 3/40, 2/25, 17/200, 171/2000, 43/500, 173/2000, 87/1000, 7/80, 11/125, 177/2000, 89/1000, 179/2000, 9/100, 19/200, 1/10, 21/200, 11/100, 23/200, 3/25, 1/8, 13/100, 27/200, 7/50, 29/200, 3/20, 31/200, 4/25, 33/200, 17/100, 7/40, 9/50, 37/200, 19/100, 39/200, 1/5, 1/4, 51/200, 13/50, 53/200, 27/100, 11/40, 7/25, 57/200, 29/100, 59/200, 3/10, 61/200, 31/100, 621/2000, 311/1000, 623/2000, 39/125, 5/16, 313/1000, 627/2000, 157/500, 629/2000, 63/200, 8/25, 13/40, 33/100, 67/200, 17/50, 69/200, 7/20, 2/5, 81/200, 41/100, 83/200, 21/50, 17/40, 43/100, 87/200, 11/25, 89/200, 9/20, 91/200, 23/50, 93/200, 47/100, 19/40, 12/25, 97/200, 49/100, 99/200, 1/2, 3/5, 61/100, 31/50, 621/1000, 311/500, 623/1000, 78/125, 5/8, 313/500, 627/1000, 157/250, 629/1000, 63/100, 631/1000, 79/125, 633/1000, 317/500, 127/200, 159/250, 637/1000, 319/500, 639/1000, 16/25, 641/1000, 321/500, 643/1000, 161/250, 129/200, 323/500, 647/1000, 81/125, 649/1000, 13/20, 33/50, 67/100, 17/25, 69/100, 7/10, 71/100, 18/25, 73/100, 37/50, 3/4, 19/25, 77/100, 39/50, 79/100, 4/5, 81/100, 41/50, 83/100, 21/25, 17/20, 43/50, 87/100, 22/25, 89/100, 891/1000, 223/250, 893/1000, 447/500, 179/200, 112/125, 897/1000, 449/500, 899/1000, 9/10, 113/125, 2261/2500, 1131/1250, 2263/2500, 566/625, 453/500, 1133/1250, 2267/2500, 567/625, 2269/2500, 227/250, 114/125, 2281/2500, 1141/1250, 2283/2500, 571/625, 457/500, 1143/1250, 2287/2500, 572/625, 2289/2500, 229/250, 23/25, 231/250, 116/125, 233/250, 117/125, 47/50, 118/125, 237/250, 119/125, 239/250, 24/25, 241/250, 121/125, 243/250, 2431/2500, 608/625, 2433/2500, 1217/1250, 487/500, 609/625, 2437/2500, 1219/1250, 2439/2500, 122/125, 49/50}, prec];


(* Internal read-only alias used by the refinement logic.  Update xSamples
   above—not this alias—when pasting the next iteration's emitted block. *)
currentSamplePoints = xSamples;

mgap = SetPrecision[166/100, prec];
nulllist = {51, -1, -1, -1};

NlistAll[z_, J_] := {(2-J (7+J))/(2 z^2),(20-J (7+J) (-13+J (7+J)))/(20 z^3),-(((-12+J (7+J)) (30+J (7+J) (-23+J (7+J))))/(360 z^4)),1/z^5-((-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J))))/(10080 z^5),-((J (7+J) (-23+J (7+J)))/(20 z^5)),1/z^6-((-4+J) (-2+J) J (7+J) (9+J) (11+J) (-34+J (5+J)) (-20+J (9+J)))/(403200 z^6),-((J (7+J) (540+J (7+J) (-53+J (7+J))))/(360 z^6)),1/z^7-((-4+J) (-2+J) J (7+J) (9+J) (11+J) (-62+J (7+J)) (-48+J (7+J)) (-15+J (7+J)))/(21772800 z^7),-((J (7+J) (-31032+J (7+J) (3024+(-7+J) J (7+J) (14+J))))/(10080 z^7)),1/z^8-((-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (-50+J (7+J)) (960+J (7+J) (-83+J (7+J))))/(1524096000 z^8),-((J (7+J) (1578240+J (7+J) (-209056+J (7+J) (8988+J (7+J) (-160+J (7+J))))))/(403200 z^8)),-(((-2+J) J (7+J) (9+J) (-53+J (7+J)))/(360 z^8)),1/z^9-1/(134120448000 z^9) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (5001600+J (7+J) (-600160+J (7+J) (19516+J (7+J) (-240+J (7+J))))),-((J (7+J) (-131466240+J (7+J) (17720496+J (7+J) (-915804+J (7+J) (21808+J (7+J) (-241+J (7+J)))))))/(21772800 z^9)),-(((-2+J) J (7+J) (9+J) (2564+J (7+J) (-108+J (7+J))))/(5040 z^9)),1/z^10-1/(14485008384000 z^10) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (5277600+J (7+J) (-651672+J (7+J) (20980+J (7+J) (-250+J (7+J))))),-(1/(1524096000 z^10))J (7+J) (11405836800+J (7+J) (-1826254080+J (7+J) (107801568+J (7+J) (-3100260+J (7+J) (46228+J (7+J) (-343+J (7+J))))))),-(((-2+J) J (7+J) (9+J) (-216000+J (7+J) (10752+J (7+J) (-182+J (7+J)))))/(201600 z^10)),1/z^11-1/(1883051089920000 z^11) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (-797572800+J (7+J) (106776000+J (7+J) (-3946556+J (7+J) (60108+J (7+J) (-405+J (7+J)))))),-(1/(134120448000 z^11))J (7+J) (-1379752704000+J (7+J) (225934848000+J (7+J) (-14770082880+J (7+J) (486509168+J (7+J) (-8812960+J (7+J) (88928+J (7+J) (-468+J (7+J)))))))),-(((-2+J) J (7+J) (9+J) (21948480+J (7+J) (-1323000+J (7+J) (28702+J (7+J) (-277+J (7+J))))))/(10886400 z^11)),-(((-2+J) J (7+J) (9+J) (35064000+J (7+J) (-1763640+J (7+J) (31942+J (7+J) (-277+J (7+J))))))/(21772800 z^11)),1/z^12-1/(289989867847680000 z^12) (-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (-833676480+J (7+J) (114378912+J (7+J) (-4230052+J (7+J) (63448+J (7+J) (-417+J (7+J)))))),-(1/(14485008384000 z^12))J (7+J) (180894269952000+J (7+J) (-33122195274240+J (7+J) (2349433112832+J (7+J) (-85849259040+J (7+J) (1788748560+J (7+J) (-22054832+J (7+J) (158952+J (7+J) (-618+J (7+J))))))))),-(1/(762048000 z^12))(-2+J) J (7+J) (9+J) (-2563473600+J (7+J) (175893840+J (7+J) (-4726576+J (7+J) (61658+J (7+J) (-395+J (7+J)))))),-(1/(1524096000 z^12))(-2+J) J (7+J) (9+J) (-5573865600+J (7+J) (318324240+J (7+J) (-6824476+J (7+J) (71108+J (7+J) (-395+J (7+J)))))),1/z^13-1/(52198176212582400000 z^13) (-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (172008748800+J (7+J) (-25008497280+J (7+J) (1017235728+J (7+J) (-17785440+J (7+J) (152128+J (7+J) (-628+J (7+J))))))),-(1/(1883051089920000 z^13))J (7+J) (-30370283513856000+J (7+J) (5685348310272000+J (7+J) (-431278994177280+J (7+J) (17117281051776+J (7+J) (-396920087856+J (7+J) (5654064144+J (7+J) (-50058872+J (7+J) (268176+J (7+J) (-795+J (7+J)))))))))),-(1/(67060224000 z^13))(-2+J) J (7+J) (9+J) (357678604800+J (7+J) (-27684313920+J (7+J) (855078608+J (7+J) (-13637120+J (7+J) (118668+J (7+J) (-538+J (7+J))))))),-(1/(134120448000 z^13))(-2+J) J (7+J) (9+J) (942637132800+J (7+J) (-59261829120+J (7+J) (1465072608+J (7+J) (-18734520+J (7+J) (134068+J (7+J) (-538+J (7+J))))))),1/z^14-((-12+J) (-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (19+J) (178803072000+J (7+J) (-26522288640+J (7+J) (1083860832+J (7+J) (-18822392+J (7+J) (158628+J (7+J) (-642+J (7+J))))))))/(10857220652217139200000 z^14),-(1/(289989867847680000 z^14))J (7+J) (5582653171070976000+J (7+J) (-1130338166370048000+J (7+J) (90722029433978880+J (7+J) (-3853417628059776+J (7+J) (97298039791872+J (7+J) (-1547871265360+J (7+J) (15899813776+J (7+J) (-105145568+J (7+J) (431816+J (7+J) (-1001+J (7+J))))))))))),-(1/(289989867847680000 z^14))J (7+J) (18904909351071744000+J (7+J) (-3159997781508403200+J (7+J) (212491113563151360+J (7+J) (-7648220048749056+J (7+J) (164888697966912+J (7+J) (-2260267109520+J (7+J) (20284193776+J (7+J) (-119680088+J (7+J) (451836+J (7+J) (-1001+J (7+J))))))))))),-(1/(7242504192000 z^14))(-2+J) J (7+J) (9+J) (-57858694502400+J (7+J) (4940899119360+J (7+J) (-172748862720+J (7+J) (3196589328+J (7+J) (-34081280+J (7+J) (211008+J (7+J) (-708+J (7+J)))))))),-(1/(14485008384000 z^14))(-2+J) J (7+J) (9+J) (-179566979328000+J (7+J) (12489984864000+J (7+J) (-350374821120+J (7+J) (5203121328+J (7+J) (-45129680+J (7+J) (234768+J (7+J) (-708+J (7+J)))))))),1/z^15-((-12+J) (-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (19+J) (-48228426086400+J (7+J) (7471309098240+J (7+J) (-327436493664+J (7+J) (6325658024+J (7+J) (-62930772+J (7+J) (336322+J (7+J) (-917+J (7+J)))))))))/(2584018515227679129600000 z^15),-(1/(52198176212582400000 z^15))J (7+J) (-1247747313202790400000+J (7+J) (257807412952510464000+J (7+J) (-21713822416662681600+J (7+J) (976289459242567680+J (7+J) (-26434504759413888+J (7+J) (459223576873344+J (7+J) (-5286173524128+J (7+J) (40717222480+J (7+J) (-207319640+J (7+J) (668976+J (7+J) (-1238+J (7+J)))))))))))),-(1/(52198176212582400000 z^15))J (7+J) (-4427589781962547200000+J (7+J) (778134090027525120000+J (7+J) (-55587311278507929600+J (7+J) (2133844925721868800+J (7+J) (-49603448066289408+J (7+J) (744185486893824+J (7+J) (-7462993855968+J (7+J) (50767274800+J (7+J) (-232960640+J (7+J) (696696+J (7+J) (-1238+J (7+J)))))))))))),-(1/(941525544960000 z^15))(-2+J) J (7+J) (9+J) (10916797731840000+J (7+J) (-1018511499110400+J (7+J) (39142949166720+J (7+J) (-814016017152+J (7+J) (10077356168+J (7+J) (-76691252+J (7+J) (353250+J (7+J) (-907+J (7+J))))))))),-(1/(1883051089920000 z^15))(-2+J) J (7+J) (9+J) (38408970023424000+J (7+J) (-2914947731980800+J (7+J) (91204340699520+J (7+J) (-1546371339552+J (7+J) (15657273368+J (7+J) (-98663852+J (7+J) (388350+J (7+J) (-907+J (7+J))))))))),1/z^16-((-14+J) (-12+J) (-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (19+J) (21+J) (-49949104896000+J (7+J) (7864546348800+J (7+J) (-346775577120+J (7+J) (6684386328+J (7+J) (-65927780+J (7+J) (347718+J (7+J) (-933+J (7+J)))))))))/(697684999111473364992000000 z^16),-((J (7+J) (305478770084167680000000+J (7+J) (-66880855243136286720000+J (7+J) (5873152091579080704000+J (7+J) (-277205314388051865600+J (7+J) (7958521963421475840+J (7+J) (-148608950149018368+J (7+J) (1873341464128704+J (7+J) (-16222546312128+J (7+J) (96556311280+J (7+J) (-387804560+J (7+J) (1003236+J (7+J) (-1508+J (7+J))))))))))))))/(10857220652217139200000 z^16)),-((J (7+J) (1188386342061144145920000+J (7+J) (-217971526056244789248000+J (7+J) (16315470556851061555200+J (7+J) (-662248807410304512000+J (7+J) (16420500781059271680+J (7+J) (-265833327608356608+J (7+J) (2921409120924864+J (7+J) (-22249607260608+J (7+J) (118056830320+J (7+J) (-431047760+J (7+J) (1040676+J (7+J) (-1508+J (7+J))))))))))))))/(10857220652217139200000 z^16)),-(1/(144994933923840000 z^16))(-2+J) J (7+J) (9+J) (-2358577574825472000+J (7+J) (237610799759116800+J (7+J) (-9941557443736320+J (7+J) (227429331438240+J (7+J) (-3164274187792+J (7+J) (28016444448+J (7+J) (-159178952+J (7+J) (563810+J (7+J) (-1137+J (7+J)))))))))),-(1/(289989867847680000 z^16))(-2+J) J (7+J) (9+J) (-9288820773872640000+J (7+J) (765214719981696000+J (7+J) (-26297097573945600+J (7+J) (496256191980480+J (7+J) (-5707230603792+J (7+J) (41913607728+J (7+J) (-200019752+J (7+J) (613860+J (7+J) (-1137+J (7+J)))))))))),1/z^17-((-14+J) (-12+J) (-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (19+J) (21+J) (17046656188416000+J (7+J) (-2775655413504000+J (7+J) (129135434676480+J (7+J) (-2693784807936+J (7+J) (29839696528+J (7+J) (-187759616+J (7+J) (673224+J (7+J) (-1280+J (7+J))))))))))/(212096239729887902957568000000 z^17),-((J (7+J) (-87853552280843169300480000+J (7+J) (19584641136608558678016000+J (7+J) (-1783834766769185877196800+J (7+J) (87808102998808260096000+J (7+J) (-2650174153472579208192+J (7+J) (52576896713946183936+J (7+J) (-714113670487896000+J (7+J) (6790569461809664+J (7+J) (-45576671260848+J (7+J) (214685735264+J (7+J) (-693742452+J (7+J) (1463280+J (7+J) (-1813+J (7+J)))))))))))))))/(2584018515227679129600000 z^17)),-((J (7+J) (-356782416111642673152000000+J (7+J) (68177518035168076333056000+J (7+J) (-5344027626418386709708800+J (7+J) (227778687728448460369920+J (7+J) (-5972738835680399729664+J (7+J) (103209702498806706432+J (7+J) (-1225240506943814592+J (7+J) (10240675349197824+J (7+J) (-60956946250928+J (7+J) (258094604768+J (7+J) (-763939124+J (7+J) (1512784+J (7+J) (-1813+J (7+J)))))))))))))))/(2584018515227679129600000 z^17)),-(1/(26099088106291200000 z^17))(-2+J) J (7+J) (9+J) (582699698380800000000+J (7+J) (-62900634396722688000+J (7+J) (2822389999565068800+J (7+J) (-69946031480532480+J (7+J) (1070137386647136+J (7+J) (-10661654391696+J (7+J) (70646061304+J (7+J) (-309730172+J (7+J) (865536+J (7+J) (-1400+J (7+J))))))))))),-(1/(52198176212582400000 z^17))(-2+J) J (7+J) (9+J) (2521949401939968000000+J (7+J) (-223540198278644736000+J (7+J) (8341871152470912000+J (7+J) (-172808387428884480+J (7+J) (2211002533791936+J (7+J) (-18416939498496+J (7+J) (102380333104+J (7+J) (-381594272+J (7+J) (934836+J (7+J) (-1400+J (7+J))))))))))),-(1/(627683696640000 z^17))(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-1098365329152000+J (7+J) (57080097100800+J (7+J) (-1240016055840+J (7+J) (14651007048+J (7+J) (-102462260+J (7+J) (427338+J (7+J) (-993+J (7+J)))))))),1/z^18-((-16+J) (-14+J) (-12+J) (-10+J) (-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (17+J) (19+J) (21+J) (23+J) (17607094156800000+J (7+J) (-2905063330752000+J (7+J) (136020135744000+J (7+J) (-2838087378720+J (7+J) (31287521168+J (7+J) (-195162000+J (7+J) (691768+J (7+J) (-1298+J (7+J))))))))))/(72112721508161887005573120000000 z^18)};

myArgs = If[Length[$ScriptCommandLine] >= 2,
  Rest[$ScriptCommandLine], {}];

If[Length[myArgs] < 2 || Length[myArgs] > 4,
  Print["USAGE: wolframscript -file refine_sampling_test17.m ",
    "<yFile> [nPts] [minWidth] <maxTotalPoints>"];
  Quit[2]
];

yFile = myArgs[[1]];

parseNumericArgument[str_String] := Quiet @ Check[
  ToExpression @ StringReplace[
    StringTrim[str],
    RegularExpression["(?<=\\d)[eE]([+-]?\\d+)"] -> "*^$1"
  ],
  $Failed
];

maxTotalPoints = parseNumericArgument[Last[myArgs]];
nPts = If[Length[myArgs] >= 3,
  parseNumericArgument[myArgs[[2]]],
  10
];
minWidth = If[Length[myArgs] >= 4,
  parseNumericArgument[myArgs[[3]]],
  1/1000000
];
outputFile = "refined_xSamples_test17.txt";

If[nPts === $Failed || !IntegerQ[nPts] || nPts < 2,
  Print["ERROR: nPts must be an integer >= 2."]; Quit[2]];
If[minWidth === $Failed || !NumericQ[minWidth] || minWidth < 0,
  Print["ERROR: minWidth must be a nonnegative number."]; Quit[2]];
If[maxTotalPoints === $Failed || !IntegerQ[maxTotalPoints] ||
   maxTotalPoints < 1,
  Print["ERROR: maxTotalPoints must be a positive integer."]; Quit[2]];

minWidth = SetPrecision[minWidth, prec];
currentSamplePoints = Sort[DeleteDuplicates[currentSamplePoints]];

If[!FileExistsQ[yFile],
  Print["ERROR: y.txt not found: ", yFile]; Quit[2]];

yRaw = Select[
  ReadList[yFile, String],
  StringLength[StringTrim[#]] > 0 &&
    !StringStartsQ[StringTrim[#], "#"] &
];
If[yRaw === {},
  Print["ERROR: y.txt is empty: ", yFile]; Quit[2]];

(* SDPB emits scientific notation with lowercase e.  Mathematica does not
   parse that notation directly, so translate digit-e-exponent to *^. *)
parseSDPBReal[str_String] := Module[{s, expr},
  s = StringReplace[
    StringTrim[str],
    RegularExpression["(?<=\\d)[eE]([+-]?\\d+)"] -> "*^$1"
  ];
  expr = Quiet[Check[ToExpression[s], $Failed]];
  If[expr === $Failed || !NumericQ[expr],
    $Failed,
    SetPrecision[N[expr, prec], prec]
  ]
];

yVec = parseSDPBReal /@ yRaw;
Do[
  If[yVec[[k]] === $Failed || !NumericQ[yVec[[k]]],
    Print["ERROR: y.txt component ", k,
      " is not numeric after parsing. Raw line: ", yRaw[[k]]];
    Quit[2]
  ],
  {k, Length[yVec]}
];

functionalCount = 2 + (nulllist[[1]] + 1);
If[Length[yVec] != functionalCount,
  Print["ERROR: y.txt has ", Length[yVec],
    " components; test17 requires ", functionalCount, "."];
  Quit[2]
];

jTiers = {
  Range[0, 1000, 2],
  Range[1500, 5000, 100],
  Range[6000, 20000, 500],
  Range[20000, 50000, 2000]
};
jValues = DeleteDuplicates[Flatten[jTiers]];

continuumVector[x_?NumericQ, j_Integer] := Module[{zv, pref},
  zv = SetPrecision[mgap + x/(1 - x), prec];
  pref = zv^18;
  N[
    Join[
      {pref*(2/zv), 0},
      pref*NlistAll[zv, j]
    ],
    prec
  ]
];

functionalValue[x_?NumericQ, j_Integer] := Module[{value},
  value = Quiet[Check[yVec . continuumVector[x, j], $Failed]];
  If[value === $Failed || !NumericQ[value],
    $Failed,
    Re[N[value, prec]]
  ]
];

worstFunctional[x_?NumericQ] := Module[
  {values, badPosition, minPosition},
  values = ParallelMap[
    functionalValue[x, #] &,
    jValues,
    Method -> "CoarsestGrained"
  ];
  badPosition = FirstPosition[values, $Failed, Missing["NotFound"]];
  If[!MissingQ[badPosition],
    Return[$Failed]
  ];
  minPosition = First[Ordering[values, 1]];
  <|"value" -> values[[minPosition]], "J" -> jValues[[minPosition]]|>
];

progressPrint[args___] := (
  Print[args];
  Scan[Flush, $Output]
);

progressPrint["=== refine_sampling_test17.m ==="];
progressPrint["  yFile             = ", yFile];
progressPrint["  nPts              = ", nPts];
progressPrint["  minWidth          = ", minWidth];
progressPrint["  maxTotalPoints    = ", maxTotalPoints];
progressPrint["  current points    = ", Length[currentSamplePoints]];
progressPrint["  unique J values   = ", Length[jValues]];
progressPrint["Launching parallel kernels..."];

LaunchKernels[];
DistributeDefinitions[
  prec, mgap, nulllist, NlistAll, yVec, jValues,
  continuumVector, functionalValue
];

progressPrint["  parallel kernels  = ", $KernelCount];
progressPrint["Starting midpoint scan of ",
  Length[currentSamplePoints] - 1, " intervals..."];

intervalData = Table[
  Module[{xa, xb, mid, width, worst},
    xa = currentSamplePoints[[i]];
    xb = currentSamplePoints[[i + 1]];
    mid = (xa + xb)/2;
    width = xb - xa;
    progressPrint["Interval ", i, "/",
      Length[currentSamplePoints] - 1, ": [", xa, ", ", xb,
      "], midpoint=", mid, " ..."];
    worst = worstFunctional[mid];
    If[worst === $Failed,
      Print["ERROR: functional evaluation failed at midpoint ", mid,
        " for interval [", xa, ", ", xb, "]."];
      Quit[2]
    ];
    progressPrint["  completed: minimum F(mid)=", worst["value"],
      " at J=", worst["J"]];
    <|
      "xa" -> xa, "xb" -> xb, "mid" -> mid, "width" -> width,
      "value" -> worst["value"], "J" -> worst["J"],
      "score" -> Abs[worst["value"]]*width
    |>
  ],
  {i, Length[currentSamplePoints] - 1}
];

narrowNegative = Select[
  intervalData,
  #["value"] < 0 && #["width"] < minWidth &
];
Do[
  Print["MIN-WIDTH SKIP: [", item["xa"], ", ", item["xb"],
    "], width=", item["width"], ", F(mid)=", item["value"],
    ", J=", item["J"]],
  {item, narrowNegative}
];

flagged = Select[
  intervalData,
  #["value"] < 0 && #["width"] >= minWidth &
];
ranked = Reverse[SortBy[flagged, #["score"] &]];

newPoints = {};
selected = {};
skipped = {};

If[flagged =!= {} && Length[currentSamplePoints] >= maxTotalPoints,
  Print["WARNING: currentSamplePoints already meets or exceeds ",
    "maxTotalPoints; refining zero intervals."];
  skipped = ranked,

  Do[
    item = ranked[[i]];
    interior = Table[
      item["xa"] + k*item["width"]/nPts,
      {k, 1, nPts - 1}
    ];
    candidate = DeleteDuplicates[
      Sort[Join[currentSamplePoints, newPoints, interior]]
    ];

    If[Length[candidate] <= maxTotalPoints,
      newPoints = Join[newPoints, interior];
      AppendTo[selected, item];
      Print["SELECTED: [", item["xa"], ", ", item["xb"],
        "], score=", item["score"], ", width=", item["width"],
        ", F(mid)=", item["value"], ", J=", item["J"],
        ", adds=", nPts - 1],

      skipped = Drop[ranked, i - 1];
      Break[]
    ],
    {i, Length[ranked]}
  ]
];

Do[
  Print["SKIPPED BY BUDGET: [", item["xa"], ", ", item["xb"],
    "], score=", item["score"], ", width=", item["width"],
    ", F(mid)=", item["value"], ", J=", item["J"]],
  {item, skipped}
];

mergedPoints = DeleteDuplicates[
  Sort[Join[currentSamplePoints, newPoints]]
];

If[Length[mergedPoints] > maxTotalPoints &&
   Length[currentSamplePoints] < maxTotalPoints,
  Print["ERROR: internal budget invariant failed: ",
    Length[mergedPoints], " > ", maxTotalPoints];
  Quit[2]
];

(* Convert the high-precision numerical grid back to exact rational numbers
   before printing.  The tolerance is far below the working accuracy, so it
   removes numerical roundoff while preserving the intended subdivision
   points.  Output therefore uses numerator/denominator form instead of
   decimal approximations. *)
rationalMergedPoints = Rationalize[
  mergedPoints,
  10^(-Floor[prec/2])
];

If[!And @@ ((IntegerQ[#] || RationalQ[#]) & /@ rationalMergedPoints),
  Print["ERROR: failed to convert every output point to an exact rational."];
  Quit[2]
];

literal = "xSamples = SetPrecision[" <>
  ToString[InputForm[rationalMergedPoints]] <>
  ", prec];";

Export[outputFile, literal <> "\n", "Text"];
Print[""];
Print["COPY THIS EXACT xSamples BLOCK INTO BOTH FILES:"];
Print[literal];
Print["Wrote literal to ", outputFile, "."];
Print["Final point count: ", Length[mergedPoints],
  " / budget ", maxTotalPoints, "."];

If[flagged === {},
  Print["STATUS: CONVERGED — zero refinable negative intervals remain."];
  Quit[0]
];

If[skipped =!= {},
  Print["STATUS: BUDGET-LIMITED, NOT CONVERGED"];
  Quit[3]
];

Print["STATUS: REFINED — all currently flagged intervals fit the budget."];
Quit[1];