function [infectCondId] = infect(gCav, Pair)
newInfectId = find(gCav.Nodes.isEmbolized);
infectCondId = zeros(1,1000);
ind = 1;
atmP = 0.1; %MPa
while ~isempty(newInfectId)
    embCond = newInfectId;
    newInfectId = zeros(1,1000);
    indnew = 1;
    for i = 1 : length(embCond)
        [eid, nid] = outedges(gCav, embCond(i));
        infectBool = gCav.Edges{eid, 'Weight'} < Pair - atmP &...
            ~ismember(nid, infectCondId);
        newInfectId(indnew : indnew + sum(infectBool) - 1) = nid(infectBool);
        indnew = indnew + length(find(infectBool));
    end
    newInfectId = unique(newInfectId(1:indnew-1));
    infectCondId(ind : ind + length(newInfectId) - 1) = newInfectId;
    ind = ind + length(newInfectId);
end
infectCondId = infectCondId(1 : ind -1);
end