function [K,Fo,gCondIt,gBipNodeIt, gCavIt] =...
    AirSeed(gCond, gBipNode, gCav,Pair,compInd,EmbNb, Ds)
%EmbNb being the number of random embolizations, in the
%beginning, for every component.
%CompInd contains the component indeces in which to embolize.
getPflag=false;
gCavIt = gCav;
gBipNodeIt = gBipNode;
gCondIt = gCond;
if EmbNb>0
    getPflag = true;
    % Identify conduits in each component.
    compBins = conncomp(gCav);
    pick = zeros(1,compInd(end)*EmbNb);
    ind = 1;
    for k = compInd
        compCond = compBins == k;
        available = find(compCond);
        i = 1;
        while i <= EmbNb
            if ~isempty(available)
                %Potential pick.
                pickPot=available(round(rand*(length(available)-1))+1);
                %Check if potential pick has already been
                %picked.
                if ismember(pickPot, pick)
                    continue
                end
                pick(ind) = pickPot;
                i = i + 1;
                ind = ind + 1;
            end
        end
    end
    [gCondIt, gBipNodeIt, gCavIt] = ...
        updateEmbolism(gCondIt, gBipNodeIt,gCavIt, pick);
end

[infectCondId] = infect(gCavIt, Pair);
if ~isempty(infectCondId)
    [gCondIt, gBipNodeIt, gCavIt] = ...
    updateEmbolism(gCondIt, gBipNodeIt,gCavIt, infectCondId);
    getPflag = true;
end
if getPflag
    %                 [~,Fo,K,~,~] = compute_hydraulics(obj,0.1,0);
    [~,Fo,K,~,~,~] =...
        compute_hydraulics(gCondIt,Pair,0.1, Ds);
else
    K=-1;Fo=-1;
end
end