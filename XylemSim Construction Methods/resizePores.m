function resizePores(obj,scheme,varargin)
validPosNum = @(x) isnumeric(x) && all(x > 0);
optInputs = inputParser;
addRequired(optInputs,'Dp',@(x) validPosNum(x) ...
    && isscalar(x) && x<200e-9 && x>1e-10);
addRequired(optInputs,'Dm',@(x) validPosNum(x) ...
    && isscalar(x) && x<200e-6 && x>1e-7);
addRequired(optInputs,'k_ASP',@(x) validPosNum(x) ...
    && isscalar(x));
addRequired(optInputs,'lam_ASP',@(x) validPosNum(x) ...
    && isscalar(x));
addRequired(optInputs,'Fc',@(x) validPosNum(x) ...
    && isscalar(x) && x<=1);
addRequired(optInputs,'Fpf',@(x) validPosNum(x) ...
    && isscalar(x) && x<=1);
addRequired(optInputs,'Fap',@(x) validPosNum(x) ...
    && isscalar(x) && x<=1);
addRequired(optInputs,'Tm',@(x) validPosNum(x) ...
    && isscalar(x));
addRequired(optInputs,'Lp',@(x) validPosNum(x) ...
    && isscalar(x));
addRequired(optInputs,'e_mean',@(x) validPosNum(x) ...
    && isscalar(x));
addRequired(optInputs,'e_cv',@(x) validPosNum(x) ...
    && isscalar(x));
parse(optInputs,varargin{:});

Dp = optInputs.Results.Dp;
Dm = optInputs.Results.Dm;
k_ASP = optInputs.Results.k_ASP;
lam_ASP = optInputs.Results.lam_ASP;
Fc = optInputs.Results.Fc;
Fpf = optInputs.Results.Fpf;
Fap = optInputs.Results.Fap;
Tm = optInputs.Results.Tm;
Lp = optInputs.Results.Lp;
e_mean = optInputs.Results.e_mean;
e_cv = optInputs.Results.e_cv;

A = k_ASP;
B = lam_ASP;
if isequal(scheme,'Deterministic') || isequal(scheme,'Stochastic')
    for i=1:length(obj.ICConnections)
        obj.ICConnections(i).A=A;
        obj.ICConnections(i).B=B;
        obj.ICConnections(i).Dp=Dp;
        obj.ICConnections(i).Dm=Dm;
        obj.ICConnections(i).Fc=Fc;
        obj.ICConnections(i).Fpf=Fpf;
        obj.ICConnections(i).Fap=Fap;
        obj.ICConnections(i).De=Dp;
        obj.ICConnections(i).Tm=Tm;
        obj.ICConnections(i).Lp=Lp;
        obj.ICConnections(i).updateMeanArea;
    end
    es = sort(abs(normrnd(e_mean,e_cv*e_mean,...
        1,length(obj.ICConnections))));
    [~,I] = sort([obj.ICConnections.Am]);
    for i=1:length(obj.ICConnections)
        obj.ICConnections(i).computeKmASP(es(I==i));
    end
    if obj.ICConnections(1).A ~= A ||...
            obj.ICConnections(1).B ~= B
        for i=1:length(obj.Conduits)
            obj.Conduits(i).addConConduitASP;
        end
    end
end
% if scheme == 2
%     ASPs = [obj.ICConnections.ASP];
%     Dpms = [obj.ICConnections.Dpm];
%     permvec = randperm(length(ASPs));
%     ASPs = ASPs(permvec);
%     Dpms = Dpms(permvec);
%     for i=1:length(obj.ICConnections)
%         obj.ICConnections(i).ASP = ASPs(i);
%         obj.ICConnections(i).Dpm = Dpms(i);
%     end
%     for i = 1:length(obj.Conduits)
%         obj.Conduits(i).addConConduits;
%     end
% end
if isequal(scheme,'Deterministic')
%     ASPs = [obj.ICConnections.ASP];
    ASPs = unique([obj.Conduits.ConConduitASP]);
    Dpms = [obj.ICConnections.Dpm];
    permvec = randperm(length(ASPs));
    ASPs = ASPs(permvec);
%     Dpms = Dpms(permvec);
    countASP = 1;
    ICCons = ICC.empty;
    countICCon = 1;
    for j = 1 : length(obj.Conduits)
        Con = obj.Conduits(j);
        
        ConConduits = Conduit.zeros(1,10);
        %         Con.ConConduitASP = zeros(1,10);
        count=1;
        for i=1:length(Con.ICConnections)
            
            if Con.ICConnections(i).ConConduits(1)==Con && ...
                    ~ismember(Con.ICConnections(i).ConConduits(2),ConConduits)
                ConConduits(count) = Con.ICConnections(i).ConConduits(2);
                if ~ismember(Con.ICConnections(i),ICCons)
                    Con.ICConnections(i).ASP = ASPs(countASP);
                    Con.ICConnections(i).Dpm = Dpms(countASP);
                    ICCons(countICCon) = Con.ICConnections(i);
                    countASP = countASP + 1;
                    countICCon = countICCon + 1;
                end
                Con.ConConduitASP(Con.ICConnections(i).ConConduits(2) == ...
                    Con.ConConduits) = Con.ICConnections(i).ASP;
                count=count+1;
            elseif Con.ICConnections(i).ConConduits(1)==Con && ...
                    ismember(Con.ICConnections(i).ConConduits(2),ConConduits)
                ind = Con.ICConnections(i).ConConduits(2) ==...
                    Con.ConConduits;
                Con.ICConnections(i).ASP = Con.ConConduitASP(ind);
            elseif Con.ICConnections(i).ConConduits(2)==Con && ...
                    ~ismember(Con.ICConnections(i).ConConduits(1),ConConduits)
                ConConduits(count) = Con.ICConnections(i).ConConduits(1);
                if ~ismember(Con.ICConnections(i),ICCons)
                    Con.ICConnections(i).ASP = ASPs(countASP);
                    Con.ICConnections(i).Dpm = Dpms(countASP);
                    ICCons(countICCon) = Con.ICConnections(i);
                    countASP = countASP + 1;
                    countICCon = countICCon + 1;
                end
                Con.ConConduitASP(Con.ICConnections(i).ConConduits(1) == ...
                    Con.ConConduits) = Con.ICConnections(i).ASP;
                count=count+1;
            elseif Con.ICConnections(i).ConConduits(2)==Con && ...
                    ismember(Con.ICConnections(i).ConConduits(1),ConConduits)
                ind = Con.ICConnections(i).ConConduits(1) == ...
                    Con.ConConduits;
                Con.ICConnections(i).ASP = Con.ConConduitASP(ind);
            end
        end
    end
end