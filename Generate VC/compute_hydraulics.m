function [Fi,Fo,Ktot,ktot,ktot_xa,ktot_sa] =...
    compute_hydraulics(sim,Pi,Po,varargin)
%Compute Flow
Pmat = sim.getPressures(Pi,Po);
[Fv,~] = sim.getFlow(Pmat);
%Calculate flow in
Fi=sum(Fv(1,:)); %m3/s
%Calculate flow out
Fo=sum(Fv(size(Fv,1)-1,:)); %m3/s

Ktot = Fo/(Pi-Po); %conductance m3/MPa.s
ktot = Ktot*(sim.Size(1)-1)*sim.Conduits(1).CEs(1).Length; %m4/MPa.s

% D1 = sim.Dcross(1);A1 = sum(pi*(D1(D1>0).^2)/4);
D2 = sim.Dcross(2,:);A2 = sum(pi*(D2(D2>0).^2)/4);
D3 = sim.Dcross(3,:);A3 = sum(pi*(D3(D3>0).^2)/4);
% D4 = sim.Dcross(4);A4 = sum(pi*(D4(D4>0).^2)/4);

ktot_xa = Ktot*(sim.Size(1)-1)*sim.Conduits(1).CEs(1).Length/...
        mean([A2,A3]);
%    mean([A1,A2,A3,A4]); %xylem area specific conductivity m2/MPa.s

if length(varargin) == 1
    VAx = varargin{1};
%     n1 = length(D1(D1>0));
    n2 = length(D2(D2>0));
    n3 = length(D3(D3>0));
%     n4 = length(D4(D4>0));
    ktot_sa = ktot*VAx/mean([n2,n3]);
%     ktot_sa = ktot*VAx/mean([n1,n2,n3,n4]);
else
    ktot_sa = NaN;
end
end