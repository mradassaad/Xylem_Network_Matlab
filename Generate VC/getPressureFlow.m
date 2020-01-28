function [Pmatrix, Fmatrix, src, sink] = getPressureFlow(gCond,Pin,Pout)
gCondLap = rmnode(gCond, find(gCond.Nodes{:,'isRedundant'} |...
    gCond.Nodes{:,'isEmbolized'}));
[wL, mWL, src, sink] = getWeightedLaplacians(gCondLap);
Pimposed = zeros(size(mWL,1),1);
srcSink = sort([src; sink]);
Pimposed(size(wL,1) + find(ismember(srcSink,src))) = Pin;
Pimposed(size(wL,1) + find(ismember(srcSink,sink))) = Pout;
PIin = mWL\Pimposed;
Pmatrix = PIin(1:size(wL,1));
%The flow is positive for (i,j) if water goes from i to j.
Fmatrix = (Pmatrix - Pmatrix' ).*adjacency(gCondLap,'weighted');
end