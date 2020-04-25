function Conduits = rmConduits(XylemNet)
            %This function removes non functional conduits, non connecting
            %clusters, and non functional ICCs
            Conduits = XylemNet.Conduits([XylemNet.Conduits.Functional]);
            FunICConnections = unique([Conduits.ICConnections]);
            
            delete(unique(...
                [obj.Conduits(~[obj.Conduits.Functional]).ICConnections]));
            delete(unique(...
                [obj.Conduits(~[obj.Conduits.Functional]).Cluster]));
            delete(obj.Conduits(~[obj.Conduits.Functional]));
            
            obj.Conduits=Conduits;
            obj.ICConnections=FunICConnections;
            obj.Clusters=obj.Clusters(isvalid(obj.Clusters));
            
            for i=1:length(obj.Conduits)
                obj.Conduits(i).addConConduits;
                obj.Conduits(i).updateConduit;
            end
            prune(obj)
end