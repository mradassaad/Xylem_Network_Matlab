function [wL, mWL, src, sink] = getWeightedLaplacians(gCond)
wA = adjacency(gCond,'weighted');
%The weighted laplacian
wL = eye(size(wA)).*sum(wA) - wA;
%network vertical distance or height
netHeigh = max(gCond.Nodes{:,'XData'});
%Identify source and sink nodes
src = find(gCond.Nodes{:,"XData"} == 1);
sink = find(gCond.Nodes{:,"XData"} == netHeigh);
srcSink = sort([src; sink]);
srcSinkMat = zeros(length(srcSink), size(wL,1));
for i = 1 : length(srcSink)
    srcSinkMat(i, srcSink(i)) = 1;
end
%Compute modified weighted Laplacian (where also the pressures at sources
%and sinks exist.
mWLTemp = [wL; srcSinkMat];
mWL = [mWLTemp [srcSinkMat'; zeros(size(srcSink,1))]];
end