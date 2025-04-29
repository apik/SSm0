SSqq = LibraryFunctionLoad["libmmaSSm0.so", "SSm0qq", {Integer, Real, Real}, Real];
SSgg = LibraryFunctionLoad["libmmaSSm0.so", "SSm0gg", {Integer, Real, Real}, Real];

beta = 0.2;
cosT = 0.3;

Table[
        Print["ep^",eo,", \\tilde{I} = ",  SSqq[eo, beta, cosT],", \\tilde{S} = ",  SSgg[eo, beta, cosT]];
       ,{eo,-3,1}]



