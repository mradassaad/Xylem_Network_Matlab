function [nodeIdx, condIdx] = redundancyCheck(gCond, gBipNode, gCav)
% This function returns the nodes and conduits that are redundant (through
% which no flow goes even though they are conducting (not isolated or
% embolized)
% The maxflow algorithm needs to run on gCond because the weights of gCav
% correspond to air-seeding pressures and not hydraulic conductances.

% comps = unique(gCav.Nodes{:,'component'});
% src = cell(1, numel(comps));
% sink = cell(1, numel(comps));
% bins = cell (1, numel(comps));
netHeight = max(gCond.Nodes{:,'XData'});

src = find(gCond.Nodes{:,'XData'} == 1 & ~gCond.Nodes{:,'isEmbolized'} &...
    ~gCond.Nodes{:,'isRedundant'});
sink = find(gCond.Nodes{:,'XData'} == netHeight &...
    ~gCond.Nodes{:,'isEmbolized'} & ~gCond.Nodes{:,'isRedundant'});

gCond = addnode(gCond,2);
gCond = addedge(gCond, numnodes(gCond) - 1, src, inf);
gCond = addedge(gCond, numnodes(gCond), sink, inf);
% for i = 1 : numel(comps)
%     compNodes = gCond.Nodes{:,'component'} == comps(i);
%     bins{i} = find(compNodes);
%     src{i} = find(gCond.Nodes{:,'XData'} == 1 & compNodes &...
%         ~gCond.Nodes{:,'isEmbolized'} & ~gCond.Nodes{:,'isRedundant'});
%     sink{i} = find(gCond.Nodes{:,'XData'} == netHeight & compNodes &...
%         ~gCond.Nodes{:,'isEmbolized'} & ~gCond.Nodes{:,'isRedundant'});
% end

[~, GF] = maxflow(gCond, numnodes(gCond) - 1, numnodes(gCond));

gCond = rmnode(gCond, numnodes(gCond)-1 : numnodes(gCond));
nodesKept = unique(GF.Edges{:,'EndNodes'}(:));
nodesKept = nodesKept(1:end-2);
condsKept = gBipNode.Edges{nodesKept, 'EndNodes'};
condsKept = condsKept(:,2) - numnodes(gCond);
% nodesKept = zeros(10000,1);
% index = 1;
% for i = 1 : numel(comps)
%     for j = 1 : length(src{i})
%        for k = 1 : length(sink{i})       
%            [~, GF] = maxflow(gCond, src{i}(j), sink{i}(k));
%            temp = GF.Edges{:,'EndNodes'}(:);
%            nodesKept(index : index + length(temp) - 1) =...
%                GF.Edges{:,'EndNodes'}(:);
%            index = index + length(temp);
%        end
%     end
% end
% 
% nodesKept = unique(nodesKept(1 : index - 1));
% condsKept = zeros(size(nodesKept));
% 
% for i = 1 : length(nodesKept)
%     condsKept(i) = neighbors(gBipNode,nodesKept(i)) - numnodes(gCond);
% end
condIdx = find(~ismember(1:numnodes(gCav),condsKept));
nodeIdx = find(~ismember(1:numnodes(gCond),nodesKept));

end