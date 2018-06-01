function [potConx] = findConx(conduitArr,Rows)
    %The first 2 rows of ConInfo will contain the nodes of the conduits and
    %the third row will contain the conduit index
    ConInfo=zeros(3,1000000);
    
    %Store Conduit Nodes
    i=1;
    for k=1:length(conduitArr)
            ConInfo(1:2,i:i+size(conduitArr(k).Nodes,1)-1) =...
                transpose(conduitArr(k).Nodes);
            ConInfo(3,i:i+size(conduitArr(k).Nodes,1)-1) = k;
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
        index = ConInfo(3,k);
        if row==1 || row==Rows
            continue
        end
        %Check if there is a horizontally adjacent node that is part of
        %another conduit. Only check nodes that are spaced by one column.
        for j=k+1:length(ConInfo)
            if abs(ConInfo(2,j) - column)>1
                break
            end
            %If Nodes are truly horizontally adjacent then designate as a
            %potential ICC
            if ConInfo(1,j)==row
                if abs(ConInfo(2,j) - column)==1
                    potConx{count}=[row column;...
                        ConInfo(1,j) ConInfo(2,j);index ConInfo(3,j)];
                    count=count+1;
                end
            end
        end
    end
    potConx=potConx(1:count-1);
end