%Modelling ring-porous wood hydraulics
%This script models the impact of having an early wood-late wood 'stepwise'
%gradient in anatomical properties on whole tissue hydraulics

%% Ring-porous
rowNb = 50; colNb = 15; depNb = 200;
Pc = [0.93 0.85]; NPc = [0.96 0.923]; Pe_rad = [0.6 0.8]; Pe_tan = [0.02 0.05];
% Average vessel length is given by 1 + Pc / (1 - Pc). 0.86 -> 7.14. A
% rowNb of 50 is enough
radDist = [ones(1, 6) ones(1, 9)*0.01]; Lce = 2.88e-3;
% Dc = 24.6e-6; Dc_cv = 0.00012;
Dc = 21.8e-6; Dc_cv = 0.12;
Dp = 32.0e-9; Dm = 6.3e-6; A = 27; B = 3.1;
fc = 0.31; fpf = 0.7;
fap = 0.06;
e_mean = 0.015;
e_cv = 0.35;
Tm = 234e-9;
Lp = 654e-9;
ASPcalcmethod = 'Pore';

rp = XylemNet(rowNb,colNb,depNb,Pc,NPc,Pe_rad,Pe_tan,radDist, Lce,Dc,Dc_cv,Dp,Dm,A,B,fc,fpf,fap,...
      e_mean,e_cv,Tm,Lp,ASPcalcmethod);
% conncomp(Agg.gCav);
[gi_val, CN, CG, VAx] = gi(rp)
viewCrossSection(rp, 20);

jt = 6; stepSize = 0.3; %MPa 
itmax = 15;  inEmb = 2; 
AvgBool = false;
[mVCRp.PLC,mVCRp.PV,mVCRp.Pressures,mVCRp.Kmax,mVCRp.Kini,mVCRp.voltot,mVCRp.embNb,...
    mVCRp.ncNb,Emb,Pen] = cavitation_process(rp,jt,stepSize,itmax,inEmb,...
    AvgBool);


% AggVC = readtable('..\Acer VCs.xlsx', 'Sheet',2);
figure
% plot(-AggVC.P,AggVC.mean)
% hold on 
plot(mVCRp.Pressures, 100 * mVCRp.PLC)

[~,~,~,ktot,ktot_xa,~] =...
    compute_hydraulics(rp.gCond,1e6,0,...
    {rp.Dcross(1,:), rp.Dcross(2,:), rp.Dcross(3,:), rp.Dcross(4,:)},...
    VAx);

disp(strcat('k_(xa) = ', num2str(ktot_xa), ' m^2/(MPa.s) or ',...
    num2str(ktot_xa*1e6/1e3), ' mg/(kPa.s.mm)'))

disp(strcat('k_(sa) = ', num2str(ktot*VAx*1e6/CN), ' m^2/(MPa.s) or ',...
    num2str(ktot*VAx*1e6/CN*1e6/1e3), ' mg/(kPa.s.mm)'))

%% Ring-porous - no change Pe
rowNb = 50; colNb = 15; depNb = 200;
Pc = [0.93 0.85]; NPc = [0.96 0.923]; Pe_rad = [0.8 0.8]; Pe_tan = [0.05 0.05];
% Average vessel length is given by 1 + Pc / (1 - Pc). 0.86 -> 7.14. A
% rowNb of 50 is enough
radDist = [ones(1, 6) ones(1, 9)*0.01]; Lce = 2.88e-3;
% Dc = 24.6e-6; Dc_cv = 0.00012;
Dc = 21.8e-6; Dc_cv = 0.12;
Dp = 32.0e-9; Dm = 6.3e-6; A = 27; B = 3.1;
fc = 0.31; fpf = 0.7;
fap = 0.06;
e_mean = 0.015;
e_cv = 0.35;
Tm = 234e-9;
Lp = 654e-9;
ASPcalcmethod = 'Pore';

rp2 = XylemNet(rowNb,colNb,depNb,Pc,NPc,Pe_rad,Pe_tan,radDist, Lce,Dc,Dc_cv,Dp,Dm,A,B,fc,fpf,fap,...
      e_mean,e_cv,Tm,Lp,ASPcalcmethod);
% conncomp(Agg.gCav);
[gi_val, CN, CG, VAx] = gi(rp2)
viewCrossSection(rp2, 20);

jt = 6; stepSize = 0.3; %MPa 
itmax = 15;  inEmb = 2; 
AvgBool = false;
[mVCRp2.PLC,mVCRp2.PV,mVCRp2.Pressures,mVCRp2.Kmax,mVCRp2.Kini,mVCRp2.voltot,mVCRp2.embNb,...
    mVCRp2.ncNb,Emb,Pen] = cavitation_process(rp2,jt,stepSize,itmax,inEmb,...
    AvgBool);


% AggVC = readtable('..\Acer VCs.xlsx', 'Sheet',2);
figure
% plot(-AggVC.P,AggVC.mean)
% hold on 
plot(mVCRp2.Pressures, 100 * mVCRp2.PLC)

[~,~,~,ktot,ktot_xa,~] =...
    compute_hydraulics(rp2.gCond,1e6,0,...
    {rp2.Dcross(1,:), rp2.Dcross(2,:), rp2.Dcross(3,:), rp2.Dcross(4,:)},...
    VAx);

disp(strcat('k_(xa) = ', num2str(ktot_xa), ' m^2/(MPa.s) or ',...
    num2str(ktot_xa*1e6/1e3), ' mg/(kPa.s.mm)'))

disp(strcat('k_(sa) = ', num2str(ktot*VAx*1e6/CN), ' m^2/(MPa.s) or ',...
    num2str(ktot*VAx*1e6/CN*1e6/1e3), ' mg/(kPa.s.mm)'))