function [infectCondId] = infect(gCav, Pair)
newInfectId = find(gCav.Nodes.isEmbolized);
infectCondId = zeros(1,1000);
ind = 1;
atmP = 0.1; %MPa
while ~isempty(newInfectId)
    embCond = newInfectId;
%     newInfectId = zeros(1,1000);
%     indnew = 1;
    %Find edges connected to embolized nodes
    embEdgeInd = arrayfun(@(x) gCav.outedges(x), embCond,...
        'UniformOutput', false);
    embEdgeInd = unique(cat(1, embEdgeInd{:}));
    %Check if concerned edges are 'occupied' or allow bubble spread under
    %current air pressure
    infectEdgeBool = gCav.Edges{embEdgeInd, 'Weight'} < Pair - atmP;
    newInfectId = gCav.Edges{embEdgeInd(infectEdgeBool),'EndNodes'};
    newInfectId = reshape(newInfectId, 1, numel(newInfectId));
    %Store only the nodes that aren't already embolized or redundant or
    %included in the previous pass over  the while loop.
    alreadyEmb = gCav.Nodes{newInfectId, 'isEmbolized'};
    alreadyRedundant = gCav.Nodes{newInfectId, 'isRedundant'};
    newInfectId = newInfectId(~(alreadyRedundant | alreadyEmb) &...
        ~ismember(newInfectId,infectCondId)');
    infectCondId(ind : ind + length(newInfectId) - 1) =...
        newInfectId;
    ind = ind + length(newInfectId);
%     for i = 1 : length(embCond)
%         [eid, nid] = outedges(gCav, embCond(i));
%         infectBool = gCav.Edges{eid, 'Weight'} < Pair - atmP &...
%             ~ismember(nid, infectCondId);
%         newInfectId(indnew : indnew + sum(infectBool) - 1) = nid(infectBool);
%         indnew = indnew + length(find(infectBool));
%     end
%     newInfectId = unique(newInfectId(1:indnew-1));
%     infectCondId(ind : ind + length(newInfectId) - 1) = newInfectId;
%     ind = ind + length(newInfectId);
% end
% infectCondId = infectCondId(1 : ind -1);
infectCondId = infectCondId(1 : ind -1);
end