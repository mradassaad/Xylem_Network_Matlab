function [ICConnections,Clusters] = generateICCs(obj,varargin)
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

%Allocate space for clusters of conduits and ICCs
Clusters = Cluster.empty();
clusCount=1;
ICClusCount = ones(1000,1);

%parameters for lognormal distribution of pit membrane diameters
% [A,B] = WblParams(ASPpit,ASPpit_cv*ASPpit);
A = k_ASP;B = lam_ASP;
ICConnections = repmat(ICC,1,length(pickedConx));
for i=1:length(pickedConx)
    %create ICC objects
    ICConnections(i) = ICC(pickedConx{i}(1:2,:),...
        [Conduits(pickedConx{i}(3,1)), Conduits(pickedConx{i}(3,2))],...
        Size,ASPcalcmethod,Dp,Dm,A,B,Fc,Fpf,Fap,Tm,Lp);
    
    %add this ICC to the corresponding conduits
    Conduits(pickedConx{i}(3,1)).ICConnections =...
        [Conduits(pickedConx{i}(3,1)).ICConnections ICConnections(i)];
    Conduits(pickedConx{i}(3,2)).ICConnections =...
        [Conduits(pickedConx{i}(3,2)).ICConnections ICConnections(i)];
    %Input Conduits into corresponding clusters
    %This case is when neither of the two conduits belond to a cluster.
    %Here, we create a new cluster.
    if isempty(Conduits(pickedConx{i}(3,1)).Cluster) &&...
            isempty(Conduits(pickedConx{i}(3,2)).Cluster)
        %Create to cluster with the two conduits
        Clusters(clusCount) = Cluster(clusCount,Size(1),...
            [Conduits(pickedConx{i}(3,1)) Conduits(pickedConx{i}(3,2))]);
        Clusters(clusCount).ICConnections = repmat(ICC,1,10000);
        %Add the ICC
        Clusters(clusCount).ICConnections(ICClusCount(clusCount)) =...
            ICConnections(i);
        ICClusCount(clusCount) = ICClusCount(clusCount)+1;
        clusCount=clusCount+1;
        
    %Next two cases is for when one of the two conduits doesn't belong in a
    %cluster, add to the other conduit's cluster
    elseif isempty(Conduits(pickedConx{i}(3,1)).Cluster)
        %tempClus is a variable that refers to the other conduit's cluster
        tempClus = Conduits(pickedConx{i}(3,2)).Cluster;
        %Add cluster-less conduit
        tempClus.AddConduits(Conduits(pickedConx{i}(3,1)));
        %Add ICC
        tempClus.ICConnections(ICClusCount(tempClus.Tag))=ICConnections(i);
        ICClusCount(tempClus.Tag) = ICClusCount(tempClus.Tag)+1;
        
    elseif isempty(Conduits(pickedConx{i}(3,2)).Cluster)
        tempClus = Conduits(pickedConx{i}(3,1)).Cluster;
        tempClus.AddConduits(Conduits(pickedConx{i}(3,2)));
        tempClus.ICConnections(ICClusCount(tempClus.Tag))=ICConnections(i);
        ICClusCount(tempClus.Tag) = ICClusCount(tempClus.Tag)+1;
    %If the two conduits belong to two different clusters, combine
    %these clusters into one
    elseif ~isequal(Conduits(pickedConx{i}(3,1)).Cluster,...
            Conduits(pickedConx{i}(3,2)).Cluster)
        
        tempClus1 = Conduits(pickedConx{i}(3,1)).Cluster;
        tempClus2 = Conduits(pickedConx{i}(3,2)).Cluster;
        %determine the number of ICCs to be added to the receiving cluster
        ICadd = ICClusCount(tempClus1.Tag)-1;
        %Copy ICCs from one cluster to the other
        tempClus2.ICConnections(ICClusCount(tempClus2.Tag):...
            ICClusCount(tempClus2.Tag)+ICadd-1) = ...
            tempClus1.ICConnections(1:ICClusCount(tempClus1.Tag)-1);
        %Copy conduits and delete one cluster
        tempClus1.IntegrateInto(tempClus2);
        ICClusCount(tempClus2.Tag) = ICClusCount(tempClus2.Tag)+ICadd;
        %Don't forget adding current ICC too
        tempClus2.ICConnections(ICClusCount(tempClus2.Tag))=ICConnections(i);
        ICClusCount(tempClus2.Tag) = ICClusCount(tempClus2.Tag)+1;
    %If the two conduits belong to the same cluster, then just add the ICC  
    elseif isequal(Conduits(pickedConx{i}(3,1)).Cluster,...
            Conduits(pickedConx{i}(3,2)).Cluster)
        
        tempClus = Conduits(pickedConx{i}(3,2)).Cluster;
        tempClus.ICConnections(ICClusCount(tempClus.Tag))=ICConnections(i);
        ICClusCount(tempClus.Tag) = ICClusCount(tempClus.Tag)+1;
    end
end

%Resize the ICC array in every cluster
for i = 1:length(Clusters)
    if isvalid(Clusters(i))
        Clusters(i).ICConnections = ...
            Clusters(i).ICConnections(1:ICClusCount(i)-1);
    end
    
end
%Keep valid clusters (those who haven't been swallowed up and deleter)
Clusters = Clusters(isvalid(Clusters));
%identify conduits in conducting clusters as functional
for i=1:length(Clusters)
    if Clusters(i).Connected
        for j=1:length(Clusters(i).Conduits)
            Clusters(i).Conduits(j).Functional=true;
        end
    end
end

% for i = 1:length(Conduits)
%     Conduits(i).addConConduits
% end
% [Conduits, ICConnections] = prune(Conduits, ICConnections);

end