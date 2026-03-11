PJ[J_, x_] := GegenbauerP[J, (D - 3)/2, x] / GegenbauerP[J, (D - 3)/2, 1];

Print[ FullSimplify[ D[ PJ[J,x], x] ] ] 