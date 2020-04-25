function [ICConnections] = generateICCs(obj,varargin)
%This function creats ICC objects, determines which conduits are connected,
%and creates clusters of connected conduits and their ICCs.
if ~isa(obj,'XylemNet')
    error('First input should be a XylemNet object')
end
p = inputParser;
validPosNum = @(x) isnumeric(x) && all(x > 0);
addRequired(p,'Dp',@(x) validPosNum(x) ...
    && isscalar(x) && x<200e-9 && x>1e-10);
addRequired(p,'Dm',@(x) validPosNum(x) ...
    && isscalar(x) && x<200e-6 && x>1e-7);
addRequired(p,'k_ASP',@(x) validPosNum(x) ...
    && isscalar(x));
addRequired(p,'lam_ASP',@(x) validPosNum(x) ...
    && isscalar(x));
addRequired(p,'Fc',@(x) validPosNum(x) ...
    && isscalar(x) && x<=1);
addRequired(p,'Fpf',@(x) validPosNum(x) ...
    && isscalar(x) && x<=1);
addRequired(p,'Fap',@(x) validPosNum(x) ...
    && isscalar(x) && x<=1);
addRequired(p,'Tm',@(x) validPosNum(x) ...
    && isscalar(x));
addRequired(p,'Lp',@(x) validPosNum(x) ...
    && isscalar(x));
addRequired(p,'ASPcalcmethod',@(x) ischar(x));
parse(p,varargin{:});

Dp = p.Results.Dp;
Dm = p.Results.Dm;
k_ASP = p.Results.k_ASP;
lam_ASP = p.Results.lam_ASP;
Fc = p.Results.Fc;
Fpf = p.Results.Fpf;
Fap = p.Results.Fap;
Tm = p.Results.Tm;
Lp = p.Results.Lp;
ASPcalcmethod = p.Results.ASPcalcmethod;

Conduits = obj.Conduits;
pickedConx = obj.pickedConx;
Size = obj.Size;

%parameters for lognormal distribution of pit membrane diameters
% [A,B] = WblParams(ASPpit,ASPpit_cv*ASPpit);
A = k_ASP;B = lam_ASP;
ICConnections = repmat(ICC,1,length(pickedConx));
for i=1:length(pickedConx)
    %create ICC objects
    ICConnections(i) = ICC(pickedConx{i}(:,1:3),...
        [Conduits(pickedConx{i}(1,4)), Conduits(pickedConx{i}(2,4))],...
        Size,ASPcalcmethod,Dp,Dm,A,B,Fc,Fpf,Fap,Tm,Lp);
    
    %add this ICC to the corresponding conduits
    Conduits(pickedConx{i}(1,4)).ICConnections =...
        [Conduits(pickedConx{i}(1,4)).ICConnections ICConnections(i)];
    Conduits(pickedConx{i}(2,4)).ICConnections =...
        [Conduits(pickedConx{i}(2,4)).ICConnections ICConnections(i)];
end

end