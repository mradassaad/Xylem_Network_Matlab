function [Conduits,SortedConduits,CBProw,CBPcol] =...
    generateConduits(NetSize, Lce, Pc, NPc, radDist)

%Pc and NPc are now arrays with two elements, where the first element
%refers to the conduits closes to the pith and the second to the conduits
%farthest radially. 
%radDist specifies how the radial profile of the vessel lengths should
%vary, radDist has the size as 'column' or NetSize(2). Where it is 1, it is
%closest to the pith, where it is 0, it is farthest radially.

%Probability of a node being the start\end of a conduit
nPcRad = (radDist*NPc(1) + (1- radDist)*NPc(2));
CondStart = arrayfun(@(x) rand(NetSize(1)+100, 1, NetSize(3))<(1-nPcRad(x)), 1: NetSize(2),...
    'UniformOutput', false);
CondStart = cat(2, CondStart{:});

PcRad = (radDist*Pc(1) + (1- radDist)*Pc(2));
CondEnd = arrayfun(@(x) rand(NetSize(1)+100, 1, NetSize(3))<(1-PcRad(x)), 1: NetSize(2),...
    'UniformOutput', false);
CondEnd = cat(2, CondEnd{:});

%Eliminate vessel number skewness toward upstream end by truncating
tempStart = zeros(1, NetSize(2), NetSize(3));
for j = 1 : NetSize(3)
    for i = 1 : NetSize(2)
        %Construct conduit at first row in this column if there exists a 1 in
        %the first 25 entries of the CondStart matrix at that column and the
        %corresponding entries of CondEnd are all 0.
        if isempty(find(CondEnd(1:50,i,j),1)) &&...
                ~isempty(find(CondStart(1:50,i,j),1))
            tempStart(1, i, j) = 1;
            %Construct conduit at first row in this column if there exists a 1 in
            %the first 25 entries of the CondStart matrix at that column that is in
            %a more advanced position than any 1 vaule in the corresponding entries
            %of CondEnd.
        elseif ~isempty(find(CondStart(1:50, i, j),1)) &&...
                find(CondEnd(1:50, i, j),1,'last') <...
                find(CondStart(1:50, i, j),1,'last')
            tempStart(1, i, j) = 1;
        end
    end
end
%Truncate CondStart and CondEnd and change their values at the initial and
%end row depending on the results of the for loop above
CondStart = CondStart(51:end-50, :, :);
CondStart(1, :, :) = tempStart;
CondStart(end, :, :) = 0;

CondEnd = CondEnd(51:end-50, :, :);
CondEnd(1, :, :) = 0;
CondEnd(end, :, :)=1;
%Clean up the contruction matrix defined as CondStart + CondEnd.*-1 where a
%conduit is initiated at a node where the value 1 is encountered or
%terminated where the value -1 is encountered.
%We need to eliminate all instances where a conduit end precedes the first
%conduit start in each column. Also, eliminate all conduit starts following
%the final conduit end in every column.
CondSE = trim(CondStart + CondEnd.*-1,1,-1,0);

%Eliminate consecutive conduit initiations or conduit terminations. In
%other words, a conduit initiation should always be followed by a conduit
%termination.
for j = 1:size(CondSE,3)
    for i = 1:size(CondSE,2)
        %Work column by column
        ttemp = CondSE(:, i, j);
        %Find positions where initiations (1) and terminations (-1) exist
        KeepR = find(ttemp);
        %Keep an initiation if it is preceded by a termination and keep a
        %termination if it is preceded by an initiation
        KeepBool = true(size(ttemp));
        KeepBool(KeepR(2:end)) = ...
            ttemp(KeepR(2:end)).*ttemp(KeepR(1:end-1)) == -1;
        %Apply changes to constructor. CondSE now only contains 0 and 1. The
        %first one is always an initiator in every column and the following 1
        %is a terminator. The 0s in between all keep on building a conduit.
        %After a termination node, all nodes with a 0 are passive and build no
        %conduits until the next initiation node is reach and the same rules
        %reapply.
        CondSE(:, i, j) = abs(CondSE(:, i, j).*KeepBool);
    end
end
%Generate conduit network blueprint (CBP: conduit blueprint)
%The blueprint is useful in case we want to create a new xylem network with
% similar conduit topology
k = find(CondSE);
[CBProw, CBPcol, CBPdep] = ind2sub(size(CondSE),k);

%Generate conduits
Conduits=repmat(Conduit(),1,10000);
SortedConduits=cell(1, NetSize(2), NetSize(3));
y=1;%counts the number of Conduits inserted
i=1;%blueprint position
%Create conductance edge table:
% 
% CondEdgeTable = table('Size',[10000 5], 'VariableTypes',...
%     {'NodeEnds' 'Type' 'Weight' 'Diameter' 'Length' 'Embolized'});

while i+1<=length(CBProw)
    %If an initiation is not followed by a terminator then remove it. This
    %case should never happen after the trimming applied above but to be
    %safe, this is here.
    if CBPcol(i) ~= CBPcol(i+1)
        CBProw(i)=[];CBPcol(i)=[];
    end
    %Initialize an array of conduit elements
    CEs=repmat(CondElem(),1,NetSize(1));
    %Construct conduit elements between conduit initiator nodes and conduit
    %terminator nodes.
    u=1;
    for j=CBProw(i):CBProw(i+1)-1
        CEs(u)=CondElem(Lce,[j CBPcol(i) CBPdep(i);j+1 CBPcol(i) CBPdep(i)],...
        NetSize);
        u=u+1;
    end
    CEs=CEs(1:u-1);
    %Create conduit object containing corresponding conduit elements,
    %taking into account the size of the network and the constant length of
    %conduit elements
    Conduits(y)=Conduit(CEs,NetSize,Lce);
    %Sort conduits per column in a cell array
    SortedConduits{1,CBPcol(i),CBPdep(i)}=...
        [SortedConduits{1,CBPcol(i),CBPdep(i)} Conduits(y)];
    y=y+1;
    %jump over the conduit terminator node to the next conduit initiator
    i=i+2;
end
Conduits=Conduits(1:y-1);
end