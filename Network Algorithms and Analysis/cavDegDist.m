cavDeg = degree(sim.gCav);
% Exclude input and output nodes
cavDeg = cavDeg(~(sim.gCav.Nodes.isInlet == 1 | sim.gCav.Nodes.isOutlet == 1));
h=histogram(cavDeg,'Normalization','probability');

countDat = h.BinCounts;
% numDat = (h.BinEdges(1:end-1)+h.BinEdges(2:end))/2;
numDat = h.BinEdges(1:end-1);
