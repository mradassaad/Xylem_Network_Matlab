function [nodeIdx, condIdx] = redundancyCheck(gCond, gBipNode, gCav, method)
% method 'maxflow' returns the nodes and conduits that are redundant (through
% which no flow goes even though they are conducting (not isolated or
% embolized).
% method 'bf' returns isolated nodes and conduits

if strcmp(method,'maxflow')
    
    % The maxflow algorithm needs to run on gCond because the weights of gCav
    % correspond to air-seeding pressures and not hydraulic conductances.

    netHeight = max(gCond.Nodes{:,'XData'});
    
    src = find(gCond.Nodes{:,'XData'} == 1 & ~gCond.Nodes{:,'isEmbolized'} &...
        ~gCond.Nodes{:,'isRedundant'});
    sink = find(gCond.Nodes{:,'XData'} == netHeight &...
        ~gCond.Nodes{:,'isEmbolized'} & ~gCond.Nodes{:,'isRedundant'});
    
    gCond = addnode(gCond,2);
    gCond = addedge(gCond, numnodes(gCond) - 1, src, inf);
    gCond = addedge(gCond, numnodes(gCond), sink, inf);

    [~, GF] = maxflow(gCond, numnodes(gCond) - 1, numnodes(gCond));
    
    gCond = rmnode(gCond, numnodes(gCond)-1 : numnodes(gCond));
    nodesKept = unique(GF.Edges{:,'EndNodes'}(:));
    nodesKept = nodesKept(1:end-2);
    condsKept = gBipNode.Edges{nodesKept, 'EndNodes'};
    condsKept = condsKept(:,2) - numnodes(gCond);

    condIdx = find(~ismember(1:numnodes(gCav),condsKept));
    nodeIdx = find(~ismember(1:numnodes(gCond),nodesKept));
    
elseif strcmp(method,'bf')
    gCavTemp = gCav;
    gCavTemp.Nodes.Index = (1:height(gCavTemp.Nodes))';
    gCavTemp = rmnode(gCavTemp, find(gCavTemp.Nodes{:,'isEmbolized'}));
    bins = conncomp(gCavTemp, 'OutputForm', 'cell');
    condIdx = zeros(1,10000);
    condPoint = 1;
    nodeIdx = zeros(1,10000);
    nodePoint = 1;
    for i = 1 : length(bins)
        if ~any(gCavTemp.Nodes{bins{i},'isInlet'}) ||...
                ~any(gCavTemp.Nodes{bins{i},'isOutlet'})
            condIdx(condPoint : condPoint + length(bins{i}) - 1) =...
                gCavTemp.Nodes{bins{i}, 'Index'};
            condPoint = condPoint + length(bins{i});
            nodeRed = arrayfun(@(x) neighbors(gBipNode, x),...
                gCavTemp.Nodes{bins{i}, 'Index'} + numnodes(gCond)...
                , 'UniformOutput', false);
            nodeRed = cat(1, nodeRed{:});
            nodeIdx(nodePoint : nodePoint + length(nodeRed) - 1) = nodeRed;
            nodePoint = nodePoint + length(nodeRed);
        end
    end
    condIdx = condIdx(1:condPoint - 1);
    nodeIdx = nodeIdx(1:nodePoint - 1);
end

end