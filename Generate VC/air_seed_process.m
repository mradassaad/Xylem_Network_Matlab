function [K,Fo]=air_seed_process(obj,Pair,ClusNb,EmbNb,SorL,RankSize,VorLPP)
%EmbNb being the number of random embolizations, in the
%beginning, for every cluster
%ClusNb determines the clusters to be embolized
getPflag=false;
if EmbNb>0
    getPflag=true;
    obj.RepairAll;
    for k=ClusNb
        if ~exist('SorL','var')
            available=obj.Clusters(k).Conduits;
        elseif strcmp(VorLPP,'V')
            Dias = [obj.Clusters(k).Conduits.Diameter];
            [~,I] = sort(Dias);
            if strcmp(SorL,'S')
                available=obj.Clusters(k).Conduits(I(1:RankSize));
            elseif strcmp(SorL,'M')
                available=obj.Clusters(k).Conduits(I(floor((end-RankSize)/2):ceil((end+RankSize)/2)));
            elseif strcmp(SorL,'L')
                available=obj.Clusters(k).Conduits(I(end-RankSize+1:end));
            end
        elseif strcmp(VorLPP,'LPP')
            Dp = [obj.Clusters(k).ICConnections.Dp];
            [~,I] = sort(Dp);
            count=0;
            i=1;
            available = Conduit.empty;
            while count < RankSize
                if strcmp(SorL,'S')
                    Con1 = obj.Clusters(k).ICConnections(I(i)).ConConduits(1);
                    Con2 = obj.Clusters(k).ICConnections(I(i)).ConConduits(2);
                elseif strcmp(SorL,'L')
                    Con1 = obj.Clusters(k).ICConnections(I(end-i+1)).ConConduits(1);
                    Con2 = obj.Clusters(k).ICConnections(I(end-i+1)).ConConduits(2);
                end
                if ~ismember(Con1,available)
                    available(count+1) = Con1;
                    count = count+1;
                end
                if ~ismember(Con2,available)
                    available(count+1) = Con2;
                    count = count+1;
                end
                i = i+1;
            end
        end
        for i=1:EmbNb
            if ~isempty(available)
                pick=round(rand*(length(available)-1))+1;
                obj.Embolize(available(pick));
                available(pick)=[];
            end
        end
    end
end
[CProp]=obj.PropagateAirSeed(Pair);

while ~isempty(CProp)
    getPflag=true;
    for i=1:length(CProp)
        obj.Embolize(CProp(i))
    end
    [CProp]=obj.PropagateAirSeed(Pair);
end
if getPflag
    [~,Fo,K,~,~] = compute_hydraulics(obj,0.1,0);
else
    K=-1;Fo=-1;
end
end