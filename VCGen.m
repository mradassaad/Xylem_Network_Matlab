net = sim;
%jt: number of VCs to be generated
jt = 5;
stepSize = 0.4; %MPa
%itmax: max number of bubble pressure increases
itmax = 16; 
%inEmb: one initial embolism in every cluster
inEmb = 1; 
%AvgBool:
%True - output averaged VC
%False - output each generation
AvgBool = false;
% 

[mVC.PLC,mVC.PV,mVC.Pressures,mVC.Kmax,mVC.Kini,mVC.voltot,mVC.embNb,...
    mVC.ncNb,Emb,Pen] = cavitation_process(sim,jt,stepSize,itmax,inEmb,...
    AvgBool);