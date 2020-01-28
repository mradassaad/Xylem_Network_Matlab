function [potConx] = findConx(conduitArr,Rows)
    %The first 3 rows of ConInfo will contain the nodes of the conduits and
    %the third row will contain the conduit index
    ConInfo=zeros(4,1000000);
    
    %Store Conduit Nodes
    i=1;
    for k=1:length(conduitArr)
            ConInfo(1:3,i:i+size(conduitArr(k).Nodes,1)-1) =...
                transpose(conduitArr(k).Nodes);
            ConInfo(4,i:i+size(conduitArr(k).Nodes,1)-1) = k;
            i=i+size(conduitArr(k).Nodes,1);
    end
    ConInfo=ConInfo(:,1:i-1);
    %Find POTENTIAL pit connections
    potConx=cell(1,100000);
    count=1;
    
    for k=1:length(ConInfo)
        %For every node that isn't at the beginning or end of the segment
        row = ConInfo(1,k);
        column = ConInfo(2,k);
        depth = ConInfo(3,k);
        index = ConInfo(4,k);
        if row==1 || row==Rows
            continue
        end
        %Check if there is a horizontally adjacent node that is part of
        %another conduit. Only check nodes that are spaced by one column.
        %The maximum number of potential connections in a 3D network is 8
        %for each node.
        for j=k+1:length(ConInfo)
            if abs(ConInfo(2,j) - column)>1 && abs(ConInfo(3,j) - depth)>1
                break
            end
            %If Nodes are truly horizontally adjacent then designate as a
            %potential ICC
            if ConInfo(1,j)==row
                if (ConInfo(2,j) - column==1 &&...
                        ConInfo(3,j) - depth<=1 &&...
                        ConInfo(3,j) - depth>=0 ) ||...
                        (ConInfo(3,j) - depth==1 &&...
                         ConInfo(2,j) - column<=1 &&...
                         ConInfo(2,j) - column>=0)
                    
                    potConx{count}=[row column depth index;...
                        ConInfo(1,j) ConInfo(2,j) ConInfo(3,j) ConInfo(4,j)];
                    count=count+1;
                end
            end
        end
    end
    potConx=potConx(1:count-1);
end