function [gCondE, gBipNodeE, gCavE] = updateEmbolism(gCond, gBipNode, gCav, condIdx)

gCavE = gCav;
gBipNodeE = gBipNode;
gCondE = gCond;

gCavE.Nodes{condIdx, 'isEmbolized'} = true;
gBipNodeE.Nodes{numnodes(gCondE) + condIdx, 'isEmbolized'} = true;
gBipNodeE.Nodes{gBipNode.Edges.EndNodes(:,1),'isEmbolized'} =...
    gBipNodeE.Nodes{gBipNode.Edges.EndNodes(:,2),'isEmbolized'};
gCondE.Nodes{:,'isEmbolized'} =...
    gBipNodeE.Nodes{1 : numnodes(gCondE), 'isEmbolized'};
gCondE.Edges{:, 'isEmbolized'} =...
    gCondE.Nodes{gCondE.Edges.EndNodes(:,1),'isEmbolized'} |...
    gCondE.Nodes{gCondE.Edges.EndNodes(:,2),'isEmbolized'};
gCondE.Edges{gCondE.Edges{:, 'isEmbolized'}, 'Weight'} = 0;
% for i = 1 : numnodes(gCondE)
%     state = gBipNodeE.Nodes{neighbors(gBipNodeE, i), 'isEmbolized'};
%     gBipNodeE.Nodes{i,'isEmbolized'} = state;
%     gCondE.Nodes{i,'isEmbolized'} = state;
%     if state
%         gCondE.Edges{outedges(gCondE, i), 'isEmbolized'} = state;
%         gCondE.Edges{outedges(gCondE, i), 'Weight'} = 0;
%     end
% end

gCavE.Nodes{:,'isRedundant'} = false;
gBipNodeE.Nodes{:,'isRedundant'} = false;
gCondE.Nodes{:,'isRedundant'} = false;
gCondE.Edges{:,'isRedundant'} = false;

[~, condIdx] = redundancyCheck(gCondE, gBipNodeE, gCavE);
[gCondE, gBipNodeE, gCavE] =...
    updateRedundancy(gCondE, gBipNodeE, gCavE, condIdx);
end