function prune(net)
% This function make sure there are no dead-end conduits in a xylem network
cons = net.Conduits;
ICCs = net.ICConnections;
NConConduits = [cons.NConConduits];
flag = true;
while flag
    flag = false;
    %determine conduits to be pruned
    tbp = cons(NConConduits <= 1);
    for i = 1 : length(tbp)
        if tbp(i).Nodes(1,1) ~= 1 && tbp(i).Nodes(end,1) ~= net.Size(1)
            conICCs = tbp(i).ICConnections;
            delete(conICCs);
            %identify connected conduit and its ICCs
            if ~isempty(tbp(i).ConConduits)
                conConID = tbp(i).ConConduits;
                %remove invalid ICC handles from connected conduit
                try
                    conConID.ICConnections =...
                        conConID.ICConnections(isvalid(conConID.ICConnections));
                catch
                    disp('1')
                end
                conConID.ICCCount = length(conConID.ICConnections);
                %delete conduit to be pruned
                delete(tbp(i));
                %reset the other fields of connected conduit
                conConID.ConConduits =...
                    conConID.ConConduits(isvalid(conConID.ConConduits));
                conConID.NConConduits = length(conConID.ConConduits);
            end
            flag = true;
        end
    end
    cons = cons(isvalid(cons));
    ICCs = ICCs(isvalid(ICCs));
%     %Also prune conduits who have only two ICCs at the same position
%     ICCCount = [cons.ICCCount];
%     %determine conduits to be potentially pruned (tbpp)
%     tbpp = cons(ICCCount == 2);
%     for i = 1:length(tbpp)
%         %It is necessary to recheck whether that conduit lost an ICC during
%         %the deletion of one before it, if so it will be pruned on the next
%         %iteration of pruning for conduits with 1 Connected conduits
%         if tbpp(i).ICCCount<2
%             continue
%         else
%             ICC1 = tbpp(i).ICConnections(1);
%             ICC2 = tbpp(i).ICConnections(2);
%         end
%         if ICC1.Nodes(1,1) == ICC2.Nodes(1,1) && tbpp(i).Nodes(1,1) ~= 1 &&...
%                 tbpp(i).Nodes(end,1) ~= net.Size(1)
%             delete([ICC1 ICC2])
%             %identify connected conduit and its ICCs
%             conConID = tbpp(i).ConConduits;
%             %delete conduit to be pruned
%             delete(tbpp(i));
%             %remove invalid ICC handles from connected conduit
%             for j = 1:length(conConID)
%                 conConID(j).ICConnections =...
%                     conConID(j).ICConnections(isvalid(conConID(j).ICConnections));
%                 conConID(j).ICCCount = length(conConID(j).ICConnections);
%                 %reset connected conduit fields
%                 conConID(j).ConConduits =...
%                     conConID(j).ConConduits(isvalid(conConID(j).ConConduits));
%                 conConID(j).NConConduits = length(conConID(j).ConConduits);
%             end
%             flag = true;
%         end
%     end
%     cons = cons(isvalid(cons));
%     ICCs = ICCs(isvalid(ICCs));
    
    NConConduits = [cons.NConConduits];
end

cons = cons(isvalid(cons));
ICCs = ICCs(isvalid(ICCs));

net.Conduits = cons;
net.ICConnections = ICCs;
for i = 1:length(net.Clusters)
    clusCons = net.Clusters(i).Conduits;
    net.Clusters(i).Conduits = clusCons(isvalid(clusCons));
    
    clusICCs = net.Clusters(i).ICConnections;
    net.Clusters(i).ICConnections = clusICCs(isvalid(clusICCs));
end

%eliminate useless conduit elements
for i = 1:length(net.Conduits)
    %Avoid conduits with only two ICCs that are on the same row. This might
    %cause the whole cluster to become isolated.
    if ~(net.Conduits(i).ICCCount == 2 &&...
            net.Conduits(i).ICConnections(1).Nodes(1,1) ==...
            net.Conduits(i).ICConnections(2).Nodes(1,1))
        net.Conduits(i).pruneCEs
    end
end

end