function [gCondU, gBipNodeU, gCavU] = updateRedundancy(gCond, gBipNode, gCav, condIdx)
% This function takes conduits indeces to list as redundant, updates gCav,
% gBipNode, and nodes and edges in gCond
gCondU = gCond;
gBipNodeU = gBipNode;
gCavU = gCav;

gCavU.Nodes{condIdx,'isRedundant'} = true;
maxNode = numnodes(gCondU);
gBipNodeU.Nodes{maxNode + 1 : end, 'isRedundant'} =...
    gCavU.Nodes{:,'isRedundant'};
gBipNodeU.Nodes{gBipNode.Edges.EndNodes(:,1),'isRedundant'} =...
    gBipNodeU.Nodes{gBipNode.Edges.EndNodes(:,2),'isRedundant'};
gCondU.Nodes{:,'isRedundant'} =...
    gBipNodeU.Nodes{1 : numnodes(gCondU), 'isRedundant'};
gCondU.Edges{:, 'isRedundant'} =...
    gCondU.Nodes{gCondU.Edges.EndNodes(:,1),'isRedundant'} |...
    gCondU.Nodes{gCondU.Edges.EndNodes(:,2),'isRedundant'};
% for i = 1 : maxNode
%     state = gBipNodeU.Nodes{neighbors(gBipNodeU, i), 'isRedundant'};
%     gBipNodeU.Nodes{i,'isRedundant'} = state;
%     gCondU.Nodes{i, 'isRedundant'} = state;
%     if state
%         gCondU.Edges{outedges(gCondU, i), 'isRedundant'} = state;
%     end
% end

end