test18Directory = If[
  StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
  DirectoryName[ExpandFileName[$InputFileName]],
  Directory[]
];
Import[FileNameJoin[{test18Directory, "..", "SDPB.m"}]];

m1 = N[3/5, 1000];
mA = N[1/1000, 1000];
J1 = 0;
J2 = 2;
mgap = N[166/100, 1000];

(* Independent functional coordinates for the three null-sum-rule families. *)
crossNullIndices = Range[0, 19];
aaNullIndices = Range[0, 19];
bbNullIndices = Range[0, 19];

nullCount = Total[
  Length /@ {crossNullIndices, aaNullIndices, bbNullIndices}
];
list0 = ConstantArray[0, nullCount];


NABAB[n_, x_, J_, m_] := {(-m^2+x)^5/x^4,(J (7+J) (m^2-x)^4)/(4 x^3),(m^2-x)^4/x^4,-((J (7+J) (m^2-x)^3)/(4 x^3)),((-1+J) J (7+J) (8+J) (m^2-x)^2)/(40 x^2),(-m^2+x)^3/x^4,(J (7+J) (m^2-x)^2)/(4 x^3),-(((-1+J) J (7+J) (8+J) (m^2-x))/(40 x^2)),((-2+J) (-1+J) J (7+J) (8+J) (9+J))/(720 x),(m^2-x)^2/x^4,(J (7+J) (-m^2+x))/(4 x^3),((-1+J) J (7+J) (8+J))/(40 x^2),-(((-2+J) (-1+J) J (7+J) (8+J) (9+J))/(720 (m^2-x) x)),((-3+J) (-2+J) (-1+J) J Gamma[11+J])/(20160 (m^2-x)^2 Gamma[7+J]),(-m^2+x)/x^4,(J (7+J))/(4 x^3),-(((-1+J) J (7+J) (8+J))/(40 (m^2-x) x^2)),((-2+J) (-1+J) J (7+J) (8+J) (9+J))/(720 (m^2-x)^2 x),-(((-3+J) (-2+J) (-1+J) J Gamma[11+J])/(20160 (m^2-x)^3 Gamma[7+J])),((-4+J) (-3+J) (-2+J) (-1+J) J x Gamma[12+J])/(806400 (m^2-x)^4 Gamma[7+J])}[[n+1]];


NAABB[n_, x_, J_, m_]:= {(J (7+J) x^(1/4) (-4 m^2+x)^(7/4) Hypergeometric2F1[1-J,8+J,5,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/(4 Sqrt[x (-4 m^2+x)]),((-4 m^2+x)^(7/4) Hypergeometric2F1[-J,7+J,4,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/x^(3/4),((-1+J) J (7+J) (8+J) (-4 m^2+x)^(3/4) Hypergeometric2F1[2-J,9+J,6,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/(40 x^(3/4)),(J (7+J) (-4 m^2+x)^(3/4) Sqrt[x (-4 m^2+x)] Hypergeometric2F1[1-J,8+J,5,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/(4 x^(7/4)),((-4 m^2+x)^(7/4) Hypergeometric2F1[-J,7+J,4,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/x^(7/4),((-2+J) (-1+J) J (7+J) (8+J) (9+J) x^(1/4) (-4 m^2+x)^(7/4) Hypergeometric2F1[3-J,10+J,7,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/(720 (x (-4 m^2+x))^(3/2)),((-1+J) J (7+J) (8+J) (-4 m^2+x)^(3/4) Hypergeometric2F1[2-J,9+J,6,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/(40 x^(7/4)),(J (7+J) (-4 m^2+x)^(7/4) Hypergeometric2F1[1-J,8+J,5,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/(4 x^(7/4) Sqrt[x (-4 m^2+x)]),((-4 m^2+x)^(7/4) Hypergeometric2F1[-J,7+J,4,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/x^(11/4),((-3+J) (-2+J) (-1+J) J Gamma[11+J] Hypergeometric2F1[4-J,11+J,8,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/(20160 x^(7/4) (-4 m^2+x)^(1/4) Gamma[7+J]),((-2+J) (-1+J) J (7+J) (8+J) (9+J) (-4 m^2+x)^(3/4) Hypergeometric2F1[3-J,10+J,7,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/(720 x^(7/4) Sqrt[x (-4 m^2+x)]),((-1+J) J (7+J) (8+J) (-4 m^2+x)^(3/4) Hypergeometric2F1[2-J,9+J,6,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/(40 x^(11/4)),(J (7+J) (-4 m^2+x)^(3/4) Sqrt[x (-4 m^2+x)] Hypergeometric2F1[1-J,8+J,5,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/(4 x^(15/4)),((-4 m^2+x)^(7/4) Hypergeometric2F1[-J,7+J,4,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/x^(15/4),((-4+J) (-3+J) (-2+J) (-1+J) J x^(1/4) (-4 m^2+x)^(7/4) Gamma[12+J] Hypergeometric2F1[5-J,12+J,9,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/(806400 (x (-4 m^2+x))^(5/2) Gamma[7+J]),((-3+J) (-2+J) (-1+J) J Gamma[11+J] Hypergeometric2F1[4-J,11+J,8,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/(20160 x^(11/4) (-4 m^2+x)^(1/4) Gamma[7+J]),((-2+J) (-1+J) J (7+J) (8+J) (9+J) (-4 m^2+x)^(7/4) Hypergeometric2F1[3-J,10+J,7,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/(720 x^(7/4) (x (-4 m^2+x))^(3/2)),((-1+J) J (7+J) (8+J) (-4 m^2+x)^(3/4) Hypergeometric2F1[2-J,9+J,6,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/(40 x^(15/4)),(J (7+J) (-4 m^2+x)^(3/4) Sqrt[x (-4 m^2+x)] Hypergeometric2F1[1-J,8+J,5,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/(4 x^(19/4)),((-4 m^2+x)^(7/4) Hypergeometric2F1[-J,7+J,4,1/2-x/(2 Sqrt[x (-4 m^2+x)])])/x^(19/4)}[[n+1]];


NBBBB[n_, z_, J_, m_] := {z-1/2 J (7+J) z,1-1/20 J (7+J) (-13+J (7+J)),-(((-12+J (7+J)) (30+J (7+J) (-23+J (7+J))))/(360 z)),1/z^2-((-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J))))/(10080 z^2),-((J (7+J) (-23+J (7+J)))/(20 z^2)),1/z^3-((-4+J) (-2+J) J (7+J) (9+J) (11+J) (-34+J (5+J)) (-20+J (9+J)))/(403200 z^3),-((J (7+J) (540+J (7+J) (-53+J (7+J))))/(360 z^3)),1/z^4-((-4+J) (-2+J) J (7+J) (9+J) (11+J) (-62+J (7+J)) (-48+J (7+J)) (-15+J (7+J)))/(21772800 z^4),-((J (7+J) (-31032+J (7+J) (3024+(-7+J) J (7+J) (14+J))))/(10080 z^4)),1/z^5-((-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (-50+J (7+J)) (960+J (7+J) (-83+J (7+J))))/(1524096000 z^5),-((J (7+J) (1578240+J (7+J) (-209056+J (7+J) (8988+J (7+J) (-160+J (7+J))))))/(403200 z^5)),-(((-2+J) J (7+J) (9+J) (-53+J (7+J)))/(360 z^5)),1/z^6-1/(134120448000 z^6)(-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (5001600+J (7+J) (-600160+J (7+J) (19516+J (7+J) (-240+J (7+J))))),-((J (7+J) (-131466240+J (7+J) (17720496+J (7+J) (-915804+J (7+J) (21808+J (7+J) (-241+J (7+J)))))))/(21772800 z^6)),-(((-2+J) J (7+J) (9+J) (2564+J (7+J) (-108+J (7+J))))/(5040 z^6)),1/z^7-1/(14485008384000 z^7)(-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (5277600+J (7+J) (-651672+J (7+J) (20980+J (7+J) (-250+J (7+J))))),-1/(1524096000 z^7)J (7+J) (11405836800+J (7+J) (-1826254080+J (7+J) (107801568+J (7+J) (-3100260+J (7+J) (46228+J (7+J) (-343+J (7+J))))))),-(((-2+J) J (7+J) (9+J) (-216000+J (7+J) (10752+J (7+J) (-182+J (7+J)))))/(201600 z^7)),1/z^8-1/(1883051089920000 z^8)(-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (-797572800+J (7+J) (106776000+J (7+J) (-3946556+J (7+J) (60108+J (7+J) (-405+J (7+J)))))),-1/(134120448000 z^8)J (7+J) (-1379752704000+J (7+J) (225934848000+J (7+J) (-14770082880+J (7+J) (486509168+J (7+J) (-8812960+J (7+J) (88928+J (7+J) (-468+J (7+J))))))))}[[n+1]];


NAAAA[n_, x_, J_, mA_] := {((-4 mA^2+x)^(3/2) (32 mA^4+2 (-1+J) (8+J) mA^2 x-(-2+J (7+J)) x^2))/(2 x^(5/2)),1/(20 x^(7/2))Sqrt[-4 mA^2+x] (-1280 mA^6+960 mA^4 x+2 (-120+(-1+J) J (7+J) (8+J)) mA^2 x^2-(-20+J (7+J) (-13+J (7+J))) x^3),1/(360 x^(9/2) Sqrt[-4 mA^2+x])(92160 mA^8-92160 mA^6 x+34560 mA^4 x^2+2 (-2880+(-2+J) (-1+J) J (7+J) (8+J) (9+J)) mA^2 x^3-(-12+J (7+J)) (30+J (7+J) (-23+J (7+J))) x^4),1/(10080 x^(11/2) (-4 mA^2+x)^(3/2))(-10321920 mA^10+12902400 mA^8 x-6451200 mA^6 x^2+1612800 mA^4 x^3+2 (-100800+(-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J)) mA^2 x^4+(10080-(-2+J) J (7+J) (9+J) (604+J (7+J) (-52+J (7+J)))) x^5),1/(180 x^(9/2) (-4 mA^2+x)^(3/2))J (7+J) (11520 mA^8-11520 mA^6 x-4 (-936+J (7+J) (-26+J (7+J))) mA^4 x^2+2 (-216+J (7+J) (-26+J (7+J))) mA^2 x^3-9 (-23+J (7+J)) x^4),1/(403200 x^(13/2) (-4 mA^2+x)^(5/2))(1651507200 mA^12-2477260800 mA^10 x+1548288000 mA^8 x^2-516096000 mA^6 x^3+96768000 mA^4 x^4+2 (-4838400+(-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J)) mA^2 x^5-(-403200+(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-34+J (5+J)) (-20+J (9+J))) x^6),1/(2520 x^(11/2) (-4 mA^2+x)^(5/2))J (7+J) (-645120 mA^10+806400 mA^8 x-403200 mA^6 x^2-2 (-54720+J (7+J) (924+J (7+J) (-56+J (7+J)))) mA^4 x^3+(-16920+J (7+J) (924+J (7+J) (-56+J (7+J)))) mA^2 x^4-7 (540+J (7+J) (-53+J (7+J))) x^5),1/(360 x^(9/2) (-4 mA^2+x)^(5/2))J (7+J) (2304 (-1+J) (8+J) mA^8+32 (-1+J) (8+J) (-90+J (7+J)) mA^6 x-24 (-1+J) (8+J) (-54+J (7+J)) mA^4 x^2+6 (-1+J) (8+J) (-42+J (7+J)) mA^2 x^3-(540+J (7+J) (-53+J (7+J))) x^4),(-356725555200 mA^14+624269721600 mA^12 x-468202291200 mA^10 x^2+195084288000 mA^8 x^3-48771072000 mA^6 x^4+7315660800 mA^4 x^5+2 (-304819200+(-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J)) mA^2 x^6-(-21772800+(-4+J) (-2+J) J (7+J) (9+J) (11+J) (-62+J (7+J)) (-48+J (7+J)) (-15+J (7+J))) x^7)/(21772800 x^(15/2) (-4 mA^2+x)^(7/2)),1/(100800 x^(13/2) (-4 mA^2+x)^(7/2))J (7+J) (103219200 mA^12-154828800 mA^10 x+96768000 mA^8 x^2-32256000 mA^6 x^3-2 (-2833920+J (7+J) (-44976+J (7+J) (3388+J (7+J) (-100+J (7+J))))) mA^4 x^4+(1+J) (6+J) (-69120+J (7+J) (4024+J (7+J) (-106+J (7+J)))) mA^2 x^5-10 (-31032+J (7+J) (3024+(-7+J) J (7+J) (14+J))) x^6),1/(10080 x^(11/2) (-4 mA^2+x)^(7/2))J (7+J) (-258048 (-1+J) (8+J) mA^10+322560 (-1+J) (8+J) mA^8 x+32 (-1+J) (8+J) (-4500+J (7+J) (-48+J (7+J))) mA^6 x^2-24 (-1+J) (8+J) (-1140+J (7+J) (-48+J (7+J))) mA^4 x^3+6 (-1+J) (8+J) (-300+J (7+J) (-48+J (7+J))) mA^2 x^4-(-31032+J (7+J) (3024+(-7+J) J (7+J) (14+J))) x^5),(99883155456000 mA^16-199766310912000 mA^14 x+174795522048000 mA^12 x^2-87397761024000 mA^10 x^3+27311800320000 mA^8 x^4-5462360064000 mA^6 x^5+682795008000 mA^4 x^6+2 (-24385536000+(-6+J) (-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J) (13+J)) mA^2 x^7-(-1524096000+(-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (-50+J (7+J)) (960+J (7+J) (-83+J (7+J)))) x^8)/(1524096000 x^(17/2) (-4 mA^2+x)^(9/2)),(J (7+J) (-44590694400 mA^14+78033715200 mA^12 x-58525286400 mA^10 x^2+24385536000 mA^8 x^3-6096384000 mA^6 x^4-4 (-240019200+J (7+J) (2888640+J (7+J) (-248256+J (7+J) (9388+J (7+J) (-160+J (7+J)))))) mA^4 x^5+2 (-49507200+J (7+J) (2888640+J (7+J) (-248256+J (7+J) (9388+J (7+J) (-160+J (7+J)))))) mA^2 x^6-27 (1578240+J (7+J) (-209056+J (7+J) (8988+J (7+J) (-160+J (7+J))))) x^7))/(10886400 x^(15/2) (-4 mA^2+x)^(9/2)),-1/(2520 x^(11/2) (-4 mA^2+x)^(9/2))(-2+J) J (7+J) (9+J) (3584 (-1+J) (8+J) mA^10+32 (-10+J) (-1+J) (8+J) (17+J) mA^8 x-32 (-1+J) (8+J) (-100+J (7+J)) mA^6 x^2+4 (-1+J) (8+J) (-230+3 J (7+J)) mA^4 x^3-2 (-1+J) (8+J) (-65+J (7+J)) mA^2 x^4+7 (-53+J (7+J)) x^5),(-35158870720512000 mA^18+79107459121152000 mA^16 x-79107459121152000 mA^14 x^2+46146017820672000 mA^12 x^3-17304756682752000 mA^10 x^4+4326189170688000 mA^8 x^5-721031528448000 mA^6 x^6+77253378048000 mA^4 x^7+2 (-2414168064000+(-7+J) (-6+J) (-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J) (13+J) (14+J)) mA^2 x^8-(-134120448000+(-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (5001600+J (7+J) (-600160+J (7+J) (19516+J (7+J) (-240+J (7+J)))))) x^9)/(134120448000 x^(19/2) (-4 mA^2+x)^(11/2)),(J (7+J) (12485394432000 mA^16-24970788864000 mA^14 x+21849440256000 mA^12 x^2-10924720128000 mA^10 x^3+3413975040000 mA^8 x^4-682795008000 mA^6 x^5-4 (-20447769600+J (7+J) (-236718720+J (7+J) (22252608+J (7+J) (-980520+J (7+J) (21868+J (7+J) (-238+J (7+J))))))) mA^4 x^6+2 (-2158617600+J (7+J) (-236718720+J (7+J) (22252608+J (7+J) (-980520+J (7+J) (21868+J (7+J) (-238+J (7+J))))))) mA^2 x^7-35 (-131466240+J (7+J) (17720496+J (7+J) (-915804+J (7+J) (21808+J (7+J) (-241+J (7+J)))))) x^8))/(762048000 x^(17/2) (-4 mA^2+x)^(11/2)),1/(50400 x^(13/2) (-4 mA^2+x)^(11/2))(-2+J) J (7+J) (9+J) (286720 (-1+J) (8+J) mA^12-430080 (-1+J) (8+J) mA^10 x-16 (-1+J) (8+J) (-15480+J (7+J) (-74+J (7+J))) mA^8 x^2+16 (-1+J) (8+J) (-4280+J (7+J) (-74+J (7+J))) mA^6 x^3-6 (-1+J) (8+J) (-1480+J (7+J) (-74+J (7+J))) mA^4 x^4+(-1+J) (8+J) (-360+J (7+J) (-74+J (7+J))) mA^2 x^5-10 (2564+J (7+J) (-108+J (7+J))) x^6),(15188632151261184000 mA^20-37971580378152960000 mA^18 x+42718027925422080000 mA^16 x^2-28478685283614720000 mA^14 x^3+12459424811581440000 mA^12 x^4-3737827443474432000 mA^10 x^5+778714050723840000 mA^8 x^6-111244864389120000 mA^6 x^7+10429206036480000 mA^4 x^8+2 (-289700167680000+(-8+J) (-7+J) (-6+J) (-5+J) (-4+J) (-3+J) (-2+J) (-1+J) J (7+J) (8+J) (9+J) (10+J) (11+J) (12+J) (13+J) (14+J) (15+J)) mA^2 x^9-(-14485008384000+(-8+J) (-6+J) (-4+J) (-2+J) J (7+J) (9+J) (11+J) (13+J) (15+J) (5277600+J (7+J) (-651672+J (7+J) (20980+J (7+J) (-250+J (7+J)))))) x^10)/(14485008384000 x^(21/2) (-4 mA^2+x)^(13/2)),(J (7+J) (-2197429420032000 mA^18+4944216195072000 mA^16 x-4944216195072000 mA^14 x^2+2884126113792000 mA^12 x^3-1081547292672000 mA^10 x^4+270386823168000 mA^8 x^5-45064470528000 mA^6 x^6-2 (-2501346355200+J (7+J) (24088008960+J (7+J) (-2417474304+J (7+J) (118343568+J (7+J) (-3123584+J (7+J) (45192+J (7+J) (-336+J (7+J)))))))) mA^4 x^7+(-388949299200+J (7+J) (24088008960+J (7+J) (-2417474304+J (7+J) (118343568+J (7+J) (-3123584+J (7+J) (45192+J (7+J) (-336+J (7+J)))))))) mA^2 x^8-22 (11405836800+J (7+J) (-1826254080+J (7+J) (107801568+J (7+J) (-3100260+J (7+J) (46228+J (7+J) (-343+J (7+J))))))) x^9))/(33530112000 x^(19/2) (-4 mA^2+x)^(13/2)),((-2+J) J (7+J) (9+J) (-123863040 (-1+J) (8+J) mA^14+216760320 (-1+J) (8+J) mA^12 x-162570240 (-1+J) (8+J) mA^10 x^2-32 (-1+J) (8+J) (-2196000+J (7+J) (5760+J (7+J) (-134+J (7+J)))) mA^8 x^3+32 (-1+J) (8+J) (-608400+J (7+J) (5760+J (7+J) (-134+J (7+J)))) mA^6 x^4-12 (-1+J) (8+J) (-290880+J (7+J) (5760+J (7+J) (-134+J (7+J)))) mA^4 x^5+2 (-1+J) (8+J) (-185040+J (7+J) (5760+J (7+J) (-134+J (7+J)))) mA^2 x^6-27 (-216000+J (7+J) (10752+J (7+J) (-182+J (7+J)))) x^7))/(5443200 x^(15/2) (-4 mA^2+x)^(13/2))}[[n+1]];


NCross[n_, z_, J_] := {
  {0, 0, -1/2 NAABB[n, z, J, mA]},
  {0, NABAB[n, z, J, mA], 0},
  {-1/2 NAABB[n, z, J, mA], 0, 0}
};

NAAMatrix[n_, z_, J_] := {
  {NAAAA[n, z, J, mA], 0, 0},
  {0, 0, 0},
  {0, 0, 0}
};

NBBMatrix[n_, z_, J_] := {
  {0, 0, 0},
  {0, 0, 0},
  {0, 0, NBBBB[n, z, J, mA]}
};

(* D=10 zero-subtracted positive moment for BB -> BB. *)
g0BBWeight[x_] :=
  2 x^2;

(* In the block-diagonal basis {AA, AB, BB}, a standard universal spin-2
   couples equally to the two neutral channels and not to AB. *)
universalSpin2Direction = {1, 0, 1};

contractUniversalSpin2[matrix_] :=
  universalSpin2Direction . matrix . universalSpin2Direction;

neutralBlock[matrix_] := matrix[[{1, 3}, {1, 3}]];


(* ---------------------------------------------------------------------- *)
(* Large-spin diagnostics (not yet imposed as PMP constraints)             *)
(* ---------------------------------------------------------------------- *)

ClearAll[
  phaseABAB10,
  aabbFixedMassGrowthBase10,
  fixedMassLargeJData10,
  largeJImpactWeights10
];

phaseABAB10[x_, m_] := (x - m^2)^7/x^4;


(* At fixed x > 4 mass^2, the AABB Gegenbauer polynomial is evaluated
   outside [-1,1] and grows exponentially with spin. *)
aabbFixedMassGrowthBase10[x_, mass_] :=
  (Sqrt[x] + 2 mass)/Sqrt[x - 4 mass^2];

fixedMassLargeJData10[x_, mass_] := <|
  "ABABCoefficientAfterJ4Scaling" -> 1/(40 x^2),
  "AABBGrowthBasePerSpin" -> aabbFixedMassGrowthBase10[x, mass],
  "AABBEvenSpinRatio" -> aabbFixedMassGrowthBase10[x, mass]^2
|>;

(* Joint large-J/large-mass limit with b = J/Sqrt[x] fixed.  If the
   physical convention is bPhysical = 2 J/Sqrt[x], use b -> bPhysical/2. *)
largeJImpactWeights10[b_, mass_] := <|
  "ABAB" -> b^4/40,
  "AABB" -> b^4 Hypergeometric0F1[6, mass^2 b^2]/40
|>;


LaunchKernels[];


PMP2SDP[datfile_, prec_:600] := Module[
    {
        npts, phiSamples, massSamples, Jmax,
        evenSpinSamples, oddSpinSamples,
        Poly, Poly1st, PolyABOdd, Poly2nd, pols, norm, obj,
        functionalCount, expectedBlocks
    },
    If[! TrueQ[Min[m1, 1, mgap] > 4 mA^2],
      Print["Invalid spectrum: every sampled pole must satisfy z > 4 mA^2."];
      Abort[]
    ];

    (* Paper eq. (D.1): Chebyshev nodes in the conformal angle.
       Here m = 1-mgap/z = Sin[phi/2]^2. *)
    npts = 200;
    phiSamples = N[Table[
      Pi/2 + Pi/2 Cos[(k + 1/2) Pi/npts],
      {k, 0, npts - 1}
    ], prec];
    massSamples = Sin[#/2]^2 & /@ phiSamples;

    If[
      Length[massSamples] =!= npts ||
        ! AllTrue[massSamples, 0 < # < 1 &],
      Print["Invalid mass-sampling grid."];
      Abort[]
    ];

    (* Neutral AA/BB states have even spin.  The independent AB spectral
       sector also admits odd spin.  Large-J constraints are not imposed yet. *)
    Jmax = 100;
    evenSpinSamples = Range[0, Jmax, 2];
    oddSpinSamples = Range[1, Jmax, 2];

    (* continuous spectrum *)
    Poly[j_, x_, y_] := Module[{g0, lambda22, polys},
      (* Block-diagonal basis {AA, AB, BB}: the neutral {AA, BB} sector is
         a 2 x 2 block, while AB is an independent 1 x 1 block. *)
      g0 = {{0, 0, 0}, {0, 0, 0}, {0, 0, g0BBWeight[x]}};

      (* Only the isolated second state contributes to lambda22. *)
      lambda22 = ConstantArray[0, {3, 3}];

      polys = Join[
        {g0, lambda22},
        NCross[#, x, j] & /@ crossNullIndices,
        NAAMatrix[#, x, j] & /@ aaNullIndices,
        NBBMatrix[#, x, j] & /@ bbNullIndices
      ];

      PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, y],
        Table[
            Table[polys[[k, row, column]], {k, Length[polys]}],
            {row, 3}, {column, 3}
        ]
       ]
    ];

    (* The first isolated state is neutral, so it only sees {AA, BB}. *)
    Poly1st[j_, x_, y_] := Module[{g0, lambda22, polys},
      g0 = {{0, 0, 0}, {0, 0, 0}, {0, 0, g0BBWeight[x]}};
      lambda22 = ConstantArray[0, {3, 3}];

      polys = neutralBlock /@ Join[
        {g0, lambda22},
        NCross[#, x, j] & /@ crossNullIndices,
        NAAMatrix[#, x, j] & /@ aaNullIndices,
        NBBMatrix[#, x, j] & /@ bbNullIndices
      ];

      PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, y],
        Table[
          Table[polys[[k, row, column]], {k, Length[polys]}],
          {row, 2}, {column, 2}
        ]
      ]
    ];

    (* Odd spins occur only in the independent AB spectral sector. *)
    PolyABOdd[j_, x_, y_] := Module[{g0, lambda22, polys},
      g0 = 0;
      lambda22 = 0;

      polys = Join[
        {g0, lambda22},
        NABAB[#, x, j, mA] & /@ crossNullIndices,
        ConstantArray[0, Length[aaNullIndices] + Length[bbNullIndices]]
      ];

      PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, y],
        {{polys}}
      ]
    ];

    (* The isolated spin-2 state has the fixed universal neutral coupling
       direction {gHAA, gHAB, gHBB} proportional to {1, 0, 1}.  Its only
       nonnegative spectral variable is the bare coupling squared lambdaH2. *)
    Poly2nd[j_, x_, y_] := Module[{g0, lambdaH2, polys},
      g0 = {{0, 0, 0}, {0, 0, 0}, {0, 0, g0BBWeight[x]}};

      (* Unit coefficient means that the second functional coordinate
         normalizes the bare universal coupling squared. *)
      lambdaH2 = 1;

      polys = Join[
        {contractUniversalSpin2[g0], lambdaH2},
        contractUniversalSpin2[NCross[#, x, j]] & /@ crossNullIndices,
        contractUniversalSpin2[NAAMatrix[#, x, j]] & /@ aaNullIndices,
        contractUniversalSpin2[NBBMatrix[#, x, j]] & /@ bbNullIndices
      ];

      PositiveMatrixWithPrefactor[
        DampedRational[1, {}, 1/E, y],
        {{polys}}
      ]
    ];
    
    pols = Flatten[{
      Flatten[ N[ ParallelTable[ Poly1st[i, m1, x], {i, J1, J1, 2}], prec] ],
      Flatten[ N[ ParallelTable[ Poly2nd[i, 1, x], {i, J2, J2, 2}], prec] ],
      Flatten[ N[ ParallelTable[ Poly[i, mgap/(1-m), x], {i, evenSpinSamples}, {m, massSamples}], prec] ],
      Flatten[ N[ ParallelTable[ PolyABOdd[i, mgap/(1-m), x], {i, oddSpinSamples}, {m, massSamples}], prec] ]
    }, 1];

    expectedBlocks = 2 + npts (
      Length[evenSpinSamples] + Length[oddSpinSamples]
    );
    If[Length[pols] =!= expectedBlocks,
      Print[
        "Unexpected block count: built ", Length[pols],
        ", expected ", expectedBlocks, "."
      ];
      Abort[]
    ];

    Print[
      "Built ", Length[pols],
      " numerical PMP blocks: neutral and AB/even J = 0, 2, ..., ", Jmax,
      "; AB/odd J = 1, 3, ..., ", Jmax - 1, "."
    ];

    (* Bound the universal isolated-state coupling squared in units of
       the BB positive moment g0BB. *)
    norm = -1 * N[Flatten[{{0, 1}, list0}], prec];
    obj = -1 * N[Flatten[{{1, 0}, list0}], prec];

    functionalCount = 2 + Length[list0];
    If[Length[norm] =!= functionalCount || Length[obj] =!= functionalCount,
      Print["Objective/normalization dimension mismatch."];
      Abort[]
    ];

    Print["size of norm = ", Length[norm]];
    Print["size of obj = ", Length[obj]];
    Print["Writing ", datfile, "..."];
    WritePmpJson[datfile, SDP[obj, norm, pols], prec, getAnalyticSampleData];

    Print["Wrote ", datfile, "."]
];

PMP2SDP[FileNameJoin[{test18Directory, "n_pmp.json"}], 1000];
