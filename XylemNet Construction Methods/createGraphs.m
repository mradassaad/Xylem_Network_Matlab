function [gCond, gBipNode, gCav] = createGraphs(XylemNet)
% gCond is the graph that dictates water translocation.
% gBip is a bipartite graph that links nodes in gCond to respective
% conduits.
% gCav shows connected conduits useful for embolism propagation
% simulations.

%Gather the conduit elements and inter-conduit connections.
ces=[XylemNet.Conduits.CEs];
iccs = [XylemNet.ICConnections];

%Retrieve the linear indices of the conduit elements and ICCs.
cesNodes = [ces.LINodes];
iccsNodes = [iccs.OINodes];

%Transform linear indices to subscipts in 3D, the last subscript
%corresponds to the CE or ICC with the largest linear index. All subscripts
%between (1,1,1) and that largest index are included in XData, YData, and
%ZData. This means that some subscripts will not be associated with any CE
%or ICC element.
[XData, YData, ZData] =...
    ind2sub(XylemNet.Size, 1: max([cesNodes(1,:) iccsNodes(1,:)...
    cesNodes(2,:) iccsNodes(2,:)]));

%First elements of node pairs - s.
sCond = [cesNodes(1,:) iccsNodes(1,:)];
%Second elements of node pairs - t.
tCond = [cesNodes(2,:) iccsNodes(2,:)];
%weights will contain the conductance of element edges or connections.
weights = zeros(size(sCond));

%ceObj is an array containing the conduit element objects. It is padded
%such that ceObj and iccObj will have the same size. iccObj is an array
%containing the inter-conduit connection objects.
ceObj = [ces repmat(CondElem(), 1, size(iccsNodes,2))];
iccObj = [repmat(ICC(), 1, size(cesNodes,2)) iccs];
edgeTypes = [repmat("CE",1,size(cesNodes,2)) repmat("ICC",1,size(iccsNodes,2))];
component = zeros(size(XData));
isRedundantEdge = false(size(sCond));
isEmbolizedEdge = false(size(sCond));
isRedundantNode = false(size(XData));
isEmbolizedNode = false(size(XData));
EdgeTableCond = table([sCond' tCond'], edgeTypes', weights', ceObj' ,iccObj',...
    isRedundantEdge', isEmbolizedEdge',...
    'VariableNames', {'EndNodes' 'Type' 'Weight' 'CEObj' 'ICCObj'...
    'isRedundant' 'isEmbolized'});
%NodeTableCond is useful to store the 3D subscripts of the nodes
NodeTableCond = table(XData', YData', ZData', component', isRedundantNode',...
    isEmbolizedNode'...
    , 'VariableNames',...
    {'XData' 'YData' 'ZData' 'component' 'isRedundant' 'isEmbolized'});

%To create the bipartite graph linking nodes to conduits, first find out
%the maximum linear index. IT IS IMPORTANT TO NOTE THAT, IN WHAT FOLLOWS,
%IT IS CRUCIAL THAT THE NODE INDICES ARE THE SAME IN GCOND AND GBIPNODE.
maxNode = max([max(sCond) max(tCond)]);
sBipNode = zeros(1,maxNode);
tBipNode = zeros(1,maxNode);
isInlet = false(1,maxNode + length(XylemNet.Conduits));
isOutlet = false(1,maxNode + length(XylemNet.Conduits));
isRedundant = false(1,maxNode + length(XylemNet.Conduits));
isEmbolized = false(1,maxNode + length(XylemNet.Conduits));
component = zeros(maxNode + length(XylemNet.Conduits),1);
index = 1;
for i = 1:length(XylemNet.Conduits)
    cond = XylemNet.Conduits(i);
    sBipNode(index:index + (cond.lastOINode - cond.firstOINode)) =...
        cond.OINodes;
    tBipNode(index:index + (cond.lastOINode - cond.firstOINode)) =...
        maxNode + i;
    if any(NodeTableCond{cond.OINodes, 'XData'} == 1)
        isInlet(maxNode + i) = true;
    end
    
    if any(NodeTableCond{cond.OINodes, 'XData'} == XylemNet.Size(1))
        isOutlet(maxNode + i) = true;
    end
    index = index + (cond.lastOINode - cond.firstOINode) + 1;
end
sBipNode = sBipNode(1:index-1);
tBipNode = tBipNode(1:index-1);
nodeTypes = [repmat("Node",1,maxNode) repmat("Conduit",1,length(XylemNet.Conduits))];
cPointer = [repmat(Conduit(), 1, maxNode), XylemNet.Conduits];
EdgeTableBip = table([sBipNode' tBipNode'],'VariableNames',...
    {'EndNodes'});
NodeTableBip = table(nodeTypes', cPointer', isInlet', isOutlet', component,...
    isRedundant', isEmbolized',...
    'VariableNames', {'Type' 'ConduitObj' 'isInlet' 'isOutlet' 'component'...
    'isRedundant' 'isEmbolized'});

gCond = graph(EdgeTableCond,NodeTableCond);
gBipNode = graph(EdgeTableBip,NodeTableBip);
%Clean up memory
clear EdgeTableCond NodeTableCond EdgeTableBip NodeTableBip...
    sBipNode tBipNode cPointer isInlet isOutlet nodeTypes...
    XData YData ZData sCond tCond edgeTypes weights CEObj ICCObj maxNode...
    cesNodes iccsNodes component

idRm = find(degree(gCond)==0);
gCond = rmnode(gCond,idRm);

idRm = find(degree(gBipNode)==0);
gBipNode = rmnode(gBipNode,idRm);

%remove clusters (components) that aren't connected to inlet or outlet:

[bins, binsizes] = conncomp(gCond);

%   1)remove those that contain less nodes than the vertical extent of the
%   segment

idRm = find(~ismember(bins, find(binsizes>=XylemNet.Size(1))));
gCond = rmnode(gCond,idRm);

%Delete the same nodes in gBipNode and then delete any lingering 'conduits'
%without any nodes attached to them.
gBipNode = rmnode(gBipNode,idRm);
idRm = find(degree(gBipNode)==0);
gBipNode = rmnode(gBipNode,idRm);

%   2)only keep those that have both an inlet and an outlet

bins = conncomp(gCond,'OutputForm','cell');
compRm = zeros(1,length(bins));
index = 1;
for i = 1:length(bins)
   if ~(any(gCond.Nodes.XData(bins{i}) == 1) &&...
           any(gCond.Nodes.XData(bins{i}) == XylemNet.Size(1)))
       
       compRm(index) = i;
       index = index + 1;
   end
end

compRm = compRm(1:index-1);
gCond = rmnode(gCond, [bins{compRm}]);
gBipNode = rmnode(gBipNode, [bins{compRm}]);
idRm = find(degree(gBipNode)==0);
gBipNode = rmnode(gBipNode,idRm);

clear idRm bins binsizes compRm

% Create graph associating ICCs to their respective conduits

condTable = gCond.Edges;
idConduits = gBipNode.Nodes.Type == "Conduit";
condArray = gBipNode.Nodes{idConduits, 'ConduitObj'};
condIsInlet = gBipNode.Nodes{idConduits, 'isInlet'};
condIsOutlet = gBipNode.Nodes{idConduits, 'isOutlet'};
condIsIsolated = gBipNode.Nodes{idConduits, 'isRedundant'};
condIsEmbolized = gBipNode.Nodes{idConduits, 'isEmbolized'};
component = zeros(length(find(idConduits)),1);
NodeTableConConduit = table(condArray , condIsInlet, condIsOutlet,...
    component, condIsIsolated, condIsEmbolized,...
    'VariableNames', {'ConduitObj' 'isInlet' 'isOutlet'...
    'component' 'isRedundant' 'isEmbolized'});
funcIdx = find(condTable.Type == "ICC");
s = zeros(1, length(funcIdx));
t = zeros(1, length(funcIdx));
for i = 1:length(funcIdx)
    endNodes = condTable{funcIdx(i), 'EndNodes'}; 
    
    cond1 = gBipNode.Nodes{neighbors(gBipNode, endNodes(1)), 'ConduitObj'};
    cond2 = gBipNode.Nodes{neighbors(gBipNode, endNodes(2)), 'ConduitObj'};
    
    cond1Idx = find(NodeTableConConduit{:, 1} == cond1, 1);
    cond2Idx = find(NodeTableConConduit{:, 1} == cond2, 1);
    
    s(i) = cond1Idx;
    t(i) = cond2Idx;
end

%Reduce multigraph to simple graph
gCav = simplify(graph(s, t, zeros(1, length(s)), NodeTableConConduit));
clear NodeTableConConduit funcIdx condTable condArray condIsInlet...
    condIsOutlet i s t endNodes cond1 cond2 cond1Idx cond2Idx

% Remove conduits with degree 1 and not an inlet or an outlet.
idRmCond = find(degree(gCav) == 1 & gCav.Nodes{:,'isInlet'} == false &...
    gCav.Nodes{:,'isOutlet'} == false);
while ~isempty(idRmCond)
   gCav = rmnode(gCav, idRmCond); 
   %Need to shit the index by the number of nodes so we actually remove
   %conduits rather than nodes. Node removal follows after this line by
   %identifying which nodes in gBipNode has degree 0. Since the node
   %indeces in gBipNode and gCond are the same, then use the same idRmNode
   %for both.
   gBipNode = rmnode(gBipNode, size(gCond.Nodes,1) + idRmCond);
   idRmNode = find(degree(gBipNode) == 0);
   gBipNode = rmnode(gBipNode, idRmNode);
   gCond = rmnode(gCond, idRmNode);
   idRmCond = find(degree(gCav) == 1 & gCav.Nodes{:,'isInlet'} == false &...
    gCav.Nodes{:,'isOutlet'} == false);
end
clear idRmCond idRmNode

%Specify remaining ICCs and CEs as functional

funcIdx = gCond.Edges.Type == "CE";
funcCE = gCond.Edges{funcIdx, 'CEObj'};
for i = 1:length(funcCE)
   funcCE(i).Functional = true; 
end

funcIdx = gCond.Edges.Type == "ICC";
funcICC = gCond.Edges{funcIdx, 'ICCObj'};
for i = 1:length(funcICC)
   funcICC(i).Functional = true; 
end
%Specify remaining conduits as functional

funcIdx = gBipNode.Nodes.Type == "Conduit";
funcCond = gBipNode.Nodes{funcIdx,'ConduitObj'};
for i = 1:length(funcCond)
   funcCond(i).Functional = true; 
end

[bins, ~] = conncomp(gCav);
gCav.Nodes{:,'component'} = bins';
maxNode = length(find(~funcIdx));
gBipNode.Nodes{funcIdx, 'component'} = bins';
for i = 1 : maxNode
    gBipNode.Nodes{i,'component'} =...
        gBipNode.Nodes{neighbors(gBipNode, i), 'component'};
end
gCond.Nodes{:, 'component'} = gBipNode.Nodes{1:maxNode,'component'};

end