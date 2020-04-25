function [Fi,Fo,Ktot,ktot,ktot_xa,ktot_sa] =...
    compute_hydraulics(gCond,Pi,Po,Ds,VAx)
% Ds is a cell array where the first element is the set of conduit diameters
% at the inlet cross-section. There are 4 cross-sections, equally spaced
% between the inlet and outlet.
%Compute Flow
[~, Fmatrix, src, sink] = getPressureFlow(gCond,Pi,Po);
%Calculate flow in
Fi = sum(Fmatrix(src,src+1),'all'); %m3/s
%Calculate flow out
Fo = sum(Fmatrix(sink-1,sink),'all'); %m3/s

netHeight = max(gCond.Nodes{:,'XData'});
ces = find(gCond.Edges.Type == "CE");
ceLength = gCond.Edges{ces(1),"CEObj"}.Length;
Ktot = Fo/(Pi-Po); %conductance m3/MPa.s
ktot = Ktot*(netHeight-1)*ceLength; %m4/MPa.s

% D1 = Ds{1} ;A1 = sum(pi*(D1(D1>0).^2)/4);
D2 = Ds{2} ;A2 = sum(pi*(D2(D2>0).^2)/4);
D3 = Ds{3} ;A3 = sum(pi*(D3(D3>0).^2)/4);
% D4 = Ds{4} ;A4 = sum(pi*(D4(D4>0).^2)/4);

ktot_xa = Ktot*(netHeight-1)*ceLength/ mean([A2,A3]);
%    mean([A1,A2,A3,A4]); %xylem area specific conductivity m2/MPa.s

if exist('VAx', 'var')
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