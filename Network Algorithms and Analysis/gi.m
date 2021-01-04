function [gi_val, gi_val_cond, condNum_cross, condGroup_cross,...
    condGroup_cross_cond, VAx] = gi(xn)
%This function calculates the grouping index (GI) of the xylem network at 4
%equidistant cross-sections. The GI is defined as the number of vessels
%divided by the number of vessel groups where, in this function, they are
%defined as components. 
% gi_cross_cond and condGroup_cross_cond only keep the ICCs that are
% available on row = yNode or yNode -1 or yNode + 1. The premise being that
% when we take a section of the stem to measure grouping, we can intersect
% any of these two intervals.

yNode = round(linspace(3, xn.Size(1) - 2, 4));

condNum_cross = zeros(1, 4);
condGroup_cross = zeros(1, 4);
condGroup_cross_cond = zeros(1, 4);
gi_cross = zeros(1,4);
gi_cross_cond = zeros(1,4);
for i = 1 : length(yNode)
%     condRm =  find(xn.gCond.Nodes.XData ~= yNode(i));
%     ggCond = rmnode(xn.gCond, condRm);
    condKeep = xn.gCond.Nodes.XData == yNode(i);
    % update gCav
    ggCavRm= xn.gBipNode.Edges(condKeep, 'EndNodes');
    ggCavRm = ggCavRm.EndNodes(:,2) - height(xn.gCond.Nodes);
    ggCavRm = setdiff(1:height(xn.gCav.Nodes), ggCavRm);
    ggCav = rmnode(xn.gCav, ggCavRm);
    
    % update gCond
    xDat = xn.gCond.Nodes.XData;
    nodeKeep = xDat == yNode(i) | xDat == yNode(i) - 1 | xDat == yNode(i) + 1;
    ggCondRm = setdiff(1:height(xn.gCond.Nodes), find(nodeKeep));
    ggCond = rmnode(xn.gCond, ggCondRm);
    condGroup_cross_cond(i) = max(conncomp(ggCond));
    
    condNum_cross(i) = height(ggCav.Nodes);
    condGroup_cross(i) = max(conncomp(ggCav));
    gi_cross(i) = condNum_cross(i)/condGroup_cross(i);
    gi_cross_cond(i) = condNum_cross(i)/condGroup_cross_cond(i);
end
avgCondArea = mean([xn.Conduits.Diameter].^2 * 1e6 ); %mm^2 - of a scare
ringArea = avgCondArea * xn.Size(2) * xn.Size(3); %mm^2

VAx = mean(condNum_cross / ringArea);
gi_val = mean(gi_cross);
gi_val_cond = mean(gi_cross_cond);
condNum_cross = mean(condNum_cross);
condGroup_cross = mean(condGroup_cross);
condGroup_cross_cond = mean(condGroup_cross_cond);
end