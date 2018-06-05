rowNb = 100;
colNb = 1000;
Pe = 0.91;
NPe = 0.87;
Pc = 0.5;
Lce = 2.88e-3;
Dc = 24.6e-6;
Dc_cv = 0.12;
Dp = 32.0e-9;
Dm = 6.3e-6;
A = 16;
B = 4;
fc = 0.31;
fpf = 0.7;
fap = 0.06;
e_mean = 0.02;
e_cv = 0.35;
Tm = 234e-9;
Lp = 654e-9;
ASPcalcmethod = 'Pore';

    sim = XylemNet(rowNb,colNb,Pe,NPe,Pc,Lce,Dc,Dc_cv,Dp,Dm,A,B,fc,fpf,fap,...
        e_mean,e_cv,Tm,Lp,ASPcalcmethod);
