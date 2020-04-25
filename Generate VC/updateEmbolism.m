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

%Check redundancy. Currently commented out because I relized that maxflow
%does not give nodes with accurate representation of of the total
%equivalent resistance. It sees the weights of edges as allowing a certain
%number of "units" per unit time. For example, it sees an edge weight of 1
%as only allowing 1 unit to pass per unit time. However, the weights in our
%case are conductances, which means that they allow any number of units to
%flow in theory, depending on the imposed potential difference,
%but a smaller weight means it's harder to do so.

%Comment out following 4 lines if using bfsearch for redundancy
% gCavE.Nodes{:,'isRedundant'} = false;
% gBipNodeE.Nodes{:,'isRedundant'} = false;
% gCondE.Nodes{:,'isRedundant'} = false;
% gCondE.Edges{:,'isRedundant'} = false;

[~, condIdx] = redundancyCheck(gCondE, gBipNodeE, gCavE, 'bf');
[gCondE, gBipNodeE, gCavE] =...
    updateRedundancy(gCondE, gBipNodeE, gCavE, condIdx);

%Check if nodes are in conducting clusters by running a breadth-first
%search



end