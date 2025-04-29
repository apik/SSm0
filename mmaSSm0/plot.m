SSqq = LibraryFunctionLoad["libmmaSSm0.so", "SSm0qq", {Integer, Real, Real}, Real];
SSgg = LibraryFunctionLoad["libmmaSSm0.so", "SSm0gg", {Integer, Real, Real}, Real];

$brdr = 10^(-5);
Print["Preparing plot for ep^1 part of (1-x)*(1-x*y)*SS_qq , safety border = ", $brdr];
pltQQ = Plot3D[(1-x)(1-x*y)*SSqq[1,x, y], {x, 0 + $brdr, 1 - $brdr}, {y, -1 + $brdr, 1 - $brdr}, PlotRange -> All];
Print["Saving results: ", Export["pltQQ.pdf", pltQQ]];


Print["Preparing plot for ep^1 part of (1-x)*(1-x*y)*SS_gg , safety border = ", $brdr];
pltGG = Plot3D[(1-x)(1-x*y)*SSgg[1,x, y], {x, 0 + $brdr, 1 - $brdr}, {y, -1 + $brdr, 1 - $brdr}, PlotRange -> All];
Print["Saving results: ", Export["pltGG.pdf", pltGG]];


Exit[]
