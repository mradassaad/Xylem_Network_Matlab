rowNb = [50 100 150]; colNb = [10 15 20]; depNb = [50 75 100];
Pc = [0.86 0.86]; NPc = [0.83 0.83]; Pe_rad = [0.9 0.9]; Pe_tan = [0.02 0.02];
Lce = 2.88e-3;
Dc = 24.6e-6; Dc_cv = 0.12;
Dp = 32.0e-9; Dm = 6.3e-6; A = 20.28; B = 3.2;
fc = 0.31; fpf = 0.7;
fap = 0.06;
e_mean = 0.015;
e_cv = 0.35;
Tm = 234e-9;
Lp = 654e-9;
ASPcalcmethod = 'Pore';
Agg = XylemNet.empty(0, 3);

for i = 1 : length(rowNb)
    Agg(i) = XylemNet(rowNb(i),colNb(i),depNb(i),Pc,NPc,Pe_rad,Pe_tan,ones(1, colNb(i)), Lce,Dc,Dc_cv,Dp,Dm,A,B,fc,fpf,fap,...
        e_mean,e_cv,Tm,Lp,ASPcalcmethod);
end

jt = 10;
stepSize = 0.3; %MPa
itmax = 15; 
inEmb = 1; 
AvgBool = false;

[PLC1,PV1,Ps1,Kmax1,Kini1,voltot1,embNb1,ncNb1,Emb1,Pen1] = cavitation_process(Agg(1),jt,stepSize,itmax,inEmb,...
        AvgBool);
[PLC2,PV2,Ps2,Kmax2,Kini2,voltot2,embNb2,ncNb2,Emb2,Pen2] = cavitation_process(Agg(2),jt,stepSize,itmax,inEmb,...
        AvgBool);
[PLC3,PV3,Ps3,Kmax3,Kini3,voltot3,embNb3,ncNb3,Emb3,Pen3] = cavitation_process(Agg(3),jt,stepSize,itmax,inEmb,...
        AvgBool);