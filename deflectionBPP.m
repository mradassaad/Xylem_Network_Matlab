function [BPP,AP] = deflectionBPP(Dm,Tm,v,E,emax,Lp,Fap)

% p10 = prctile(data,10);
% p90 = prctile(data,90);

Dap = Dm*(Fap)^(1/2); %m - diameter of pit aperture
%AP - Aspiraton Pressure
AP = (1/3)*E*Tm^3*Lp/((1-v^2)*(Dm^2 - Dap^2)^2);

%BPP - Buble Propagation Pressure
BPP = emax*((3/8)*(1-v^2)*(Dm/2)^2/(E*Tm^2))^(-1); %MPa
if BPP > AP
    %eAsp is the maximum strain at aspiration
    eAsp = ((3/8)*(1-v^2)*(Dm/2)^2/(E*Tm^2))*AP; %MPa
    %eRest is the extra maximum strain the pit membrane can hold before BPP
    eRest = emax - eAsp;
    %PAsp is the addition pressure held by the membrane after aspiration
    PAsp = eRest*((3/8)*(1-v^2)*(Dap/2)^2/(E*Tm^2))^(-1); %MPa
    
    BPP = AP + PAsp;
end

end