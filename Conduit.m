classdef Conduit < handle
    properties
        CEs %Conduit Elements
        ICConnections %Inter-Conduit Connectections
        Nodes
        OINodes
        firstOINode
        lastOINode
        Functional %Logical that specifies if the conduit is part of a
        %functional cluster
        Cluster %Specifies which cluster this conduit belongs to
        Diameter %meters
        ConConduits %conduits connected to this one
        ConConduitASP %minimum ASP with each conduit
        ConConduitNpit
        ConConduitNpore
        Kce %Conductance of constituting conduit elements in %m^3/MPa.s
        ICCCount
        NConConduits
        netSize %Size of the network where this conduit exists
        isInlet
        isOutlet
    end
    
    properties (SetAccess = private)
        CECount
        Length
        Lce
    end
    
    properties (SetAccess = private)
        Embolized %Specifies whether this conduit is embolized
        NonConducting %Specifies whether a local embolism or embolism elsewhere isolates this conduit
    end
    
    properties (Dependent)
        Ac %surface area of contuit
        %ICConnections properties
        Nm %total number of pit connections in the conduit
    end
    
    methods
        
        function obj = Conduit(cellElem,netSize,Lce)
            %Check if this conduit has nodes in it
            if nargin<=0
                obj.Diameter = 0;
                obj.Length = 0;
                obj.CEs = CondElem.empty;
                obj.ICConnections = ICC.empty;
                obj.Functional = false;
                obj.Cluster = Cluster.empty();
                obj.Embolized = false;
                obj.NonConducting = false;
                obj.Kce = 0;
                obj.Lce = 0;
            else
                obj.CECount = length(cellElem);

                %Create conduit elements matrix
                obj.CEs=cellElem;
                
                %Include nodes
                %Create nodes matrix
                obj.Nodes=zeros(length(cellElem)+1,2);
                
                for k=1:length(obj.CEs)
                    obj.Nodes(k,:)=cellElem(k).Nodes(1,:);
                end
                obj.Nodes(end,:)=cellElem(end).Nodes(2,:);
                %Determine if this conduit is a water inlet for the segment
                if obj.Nodes(1,1) == 1
                    obj.isInlet = true;
                else
                    obj.isInlet = false;
                end
                %Determine if this conduit is a water outlet for the segment
                if obj.Nodes(end,1) == netSize(1)
                    obj.isOutlet = true;
                else
                    obj.isOutlet = false;
                end
                    
                %make sure all elements have the same diameter
                dia = obj.CEs(1).Diameter; %m
                for k=2:length(cellElem)
                    if obj.CEs(k).Diameter ~= dia
                        error('Elements should have the same diameter')
                    end
                end
                
                obj.Diameter = dia; %m
                obj.Length = cellElem(1).Length*obj.CECount; %m
                obj.CEs = cellElem;
                obj.Lce = Lce;
                
                obj.ICConnections=ICC.empty;
                obj.Functional=false;
                obj.Cluster=Cluster.empty();
                obj.Embolized=false;
                obj.NonConducting=false;
                obj.Kce=obj.CEs(1).Kce;
                obj.OINodes=obj.OneIndexNodes(netSize);
                obj.firstOINode=obj.OINodes(1);
                obj.lastOINode=obj.OINodes(end);
                obj.netSize = netSize;
            end
        end
        
        function updateConduit(obj)
            obj.ICCCount=length(obj.ICConnections);
            obj.NConConduits=length(obj.ConConduits);
        end
        function updateDiameter(obj,dia)
            if obj.CEs(1).Diameter == 0
                arrayfun(@(x) x.updateDiameter(dia),obj.CEs);
            end
            obj.Diameter = dia; %m
            obj.Kce=obj.CEs(1).Kce;
        end
        
        function pruneCEs(obj) 
            %Because ICConnections of every conduit are sorted according
            %to their vertical position, we can eliminate all CEs below
            %the first ICC and all CEs above the last ICC. However, check
            %that first
            
            ICCs = obj.ICConnections;
            %ICCind in this situation is the vertical (row) index of ICC
            ICCind = zeros(size(ICCs));
            for i = 1:length(ICCs)
               ICCind(i) = ICCs(i).Nodes(1); 
            end
            [~,I] = sort(ICCind);
            obj.ICConnections = obj.ICConnections(I);
            try
            firstICCInd = obj.ICConnections(1).Nodes(1,1);
            catch
                disp('1')
            end
            lastICCInd = obj.ICConnections(end).Nodes(1,1);
            %Now identify CEs to prune.
            %Find first and last vertical indices of conduit
            firstConInd = obj.Nodes(1,1);
            lastConInd = obj.Nodes(end,1);
            %We can also use the fact that all CEs are sorted according to
            %their vertical index
            if lastConInd == obj.netSize(1) && firstConInd ~= 1
                delete(obj.CEs(1:firstICCInd-firstConInd))
            elseif lastConInd ~= obj.netSize(1) && firstConInd == 1
                delete(obj.CEs(end-(lastConInd-lastICCInd-1):end))
            else
                delete(obj.CEs(1:firstICCInd-firstConInd))
                delete(obj.CEs(end-(lastConInd-lastICCInd-1):end))
            end
                %update CEs array
                obj.CEs = obj.CEs(isvalid(obj.CEs));
                
                obj.CECount = length(obj.CEs);
                
                obj.Nodes=zeros(obj.CECount+1,2);
                
                for k=1:obj.CECount
                    obj.Nodes(k,:)=obj.CEs(k).Nodes(1,:);
                end
                try
                obj.Nodes(end,:)=obj.CEs(end).Nodes(2,:);
                catch
                    disp('1')
                end
                
                obj.Length = obj.CEs(1).Length*obj.CECount; %m
                
                obj.OINodes=obj.OneIndexNodes(obj.netSize);
                obj.firstOINode=obj.OINodes(1);
                obj.lastOINode=obj.OINodes(end);
        end
        function Area = get.Ac(obj)
            Area = obj.Length*pi*obj.Diameter; %m^2
        end
        
        function number = get.Nm(obj)
            number = length(obj.ICConnections);
        end
        
        function OINodes = OneIndexNodes(obj,dimensions)
            %Convert conduit and pit matrices to equivalent positions in the new
            %numbering system
            OINodes=sub2ind(dimensions,obj.Nodes(:,1),obj.Nodes(:,2));
        end
        
        function OINodesShift = ShiftedOneIndexNodes(obj,dimensions,clusterColumns)
            %Convert conduit and pit matrices to equivalent positions in the new
            %numbering system inside clusters
            OINodesShift=dimensions(1)*((obj.Nodes(:,2)-clusterColumns(1)+1)-1)+obj.Nodes(:,1);
        end
        
%         function Input = isinput(obj)
%             %Returns 1 if this conduit is an input
%             if obj.Nodes(1,1)==1
%                 Input=1;
%             else
%                 Input=0;
%             end
%         end
%         
%         function Output = isoutput(obj,dim)
%             if obj.Nodes(end,1)==dim(1)
%                 Output=1;
%             else
%                 Output=0;
%             end
%         end
        
        function Embolize(obj)
            obj.Embolized=true;
            obj.NonConducting=true;
            [obj.ICConnections.Embolized] = deal(true);
            [obj.ICConnections.NonConducting] = deal(true);
%             arrayfun(@(x) x.Embolize,[obj.ICConnections]);
        end
        
        function Isolate(obj)
            obj.NonConducting=true;
            [obj.ICConnections.NonConducting] = deal(true);
%             arrayfun(@(x) x.Isolate,[obj.ICConnections]);
        end
        
        function Repair(obj)
            obj.Embolized=false;
            [obj.ICConnections.Embolized] = deal(false);
%             arrayfun(@(x) x.Repair,[obj.ICConnections]);
        end
        
        function Integrate(obj)
            obj.NonConducting=false;
            [obj.ICConnections.NonConducting] = deal(false);
%             arrayfun(@(x) x.Integrate,[obj.ICConnections]);
        end
        
        function addConConduits(obj)
            if isempty(obj.ICConnections)
                error('No Inter-Conduit Connections')
            end
            obj.ConConduits=Conduit.zeros(1,10);
%             obj.ConConduitASP = zeros(1,10);
            count=1;
            if ~isempty(obj.ICConnections)
                for i=1:length(obj.ICConnections)
                    if obj.ICConnections(i).ConConduits(1)==obj && ...
                            ~ismember(obj.ICConnections(i).ConConduits(2),obj.ConConduits)
                        obj.ConConduits(count) = obj.ICConnections(i).ConConduits(2);
                        %                         obj.ConConduitASP(count) = obj.ICConnections(i).ASP;
                        count=count+1;
                    elseif obj.ICConnections(i).ConConduits(1)==obj && ...
                            ismember(obj.ICConnections(i).ConConduits(2),obj.ConConduits)
%                         ind = find(obj.ICConnections(i).ConConduits(2) ==...
%                             obj.ConConduits);
                        %                         obj.ConConduitASP(ind) =...
                        %                             min([obj.ConConduitASP(ind) obj.ICConnections(i).ASP]);
                    elseif obj.ICConnections(i).ConConduits(2)==obj && ...
                            ~ismember(obj.ICConnections(i).ConConduits(1),obj.ConConduits)
                        obj.ConConduits(count) = obj.ICConnections(i).ConConduits(1);
                        %                         obj.ConConduitASP(count) = obj.ICConnections(i).ASP;
                        count=count+1;
                    elseif obj.ICConnections(i).ConConduits(2)==obj && ...
                            ismember(obj.ICConnections(i).ConConduits(1),obj.ConConduits)
%                         ind = find(obj.ICConnections(i).ConConduits(1) == ...
%                             obj.ConConduits);
                        %                         obj.ConConduitASP(ind) =...
                        %                             min([obj.ConConduitASP(ind) obj.ICConnections(i).ASP]);
                    end
                end
            end
            obj.ConConduits = obj.ConConduits(1:count-1);
%             obj.ConConduitASP = obj.ConConduitASP(1:count-1);
        end
        
        function addConConduitASP(obj)
            if isempty(obj.ICConnections)
                error('No Inter-Conduit Connections')
            end
            obj.ConConduits=Conduit.zeros(1,10);
            obj.ConConduitASP = zeros(1,10);
            obj.ConConduitNpit = zeros(1,10);
            obj.ConConduitNpore = zeros(1,10);
            count=1;
            if ~isempty(obj.ICConnections)
                for i=1:length(obj.ICConnections)
                    if obj.ICConnections(i).ConConduits(1)==obj && ...
                            ~ismember(obj.ICConnections(i).ConConduits(2),obj.ConConduits)
                        obj.ConConduits(count) = obj.ICConnections(i).ConConduits(2);
                        obj.ConConduitASP(count) = obj.ICConnections(i).ASP;
                        obj.ConConduitNpit(count) = obj.ICConnections(i).Npit;
                        obj.ConConduitNpore(count) = obj.ICConnections(i).Npore;
                        count=count+1;
                    elseif obj.ICConnections(i).ConConduits(1)==obj && ...
                            ismember(obj.ICConnections(i).ConConduits(2),obj.ConConduits)
                        ind = find(obj.ICConnections(i).ConConduits(2) ==...
                            obj.ConConduits);
                        obj.ConConduitASP(ind) =...
                            min([obj.ConConduitASP(ind) obj.ICConnections(i).ASP]);
                        obj.ConConduitNpit(ind) = obj.ConConduitNpit(ind)+...
                            obj.ICConnections(i).Npit;
                        obj.ConConduitNpore(ind) = obj.ConConduitNpore(ind)+...
                            obj.ICConnections(i).Npore;
                    elseif obj.ICConnections(i).ConConduits(2)==obj && ...
                            ~ismember(obj.ICConnections(i).ConConduits(1),obj.ConConduits)
                        obj.ConConduits(count) = obj.ICConnections(i).ConConduits(1);
                        obj.ConConduitASP(count) = obj.ICConnections(i).ASP;
                        obj.ConConduitNpit(count) = obj.ICConnections(i).Npit;
                        obj.ConConduitNpore(count) = obj.ICConnections(i).Npore;
                        count=count+1;
                    elseif obj.ICConnections(i).ConConduits(2)==obj && ...
                            ismember(obj.ICConnections(i).ConConduits(1),obj.ConConduits)
                        ind = find(obj.ICConnections(i).ConConduits(1) == ...
                            obj.ConConduits);
                        obj.ConConduitASP(ind) =...
                            min([obj.ConConduitASP(ind) obj.ICConnections(i).ASP]);
                        obj.ConConduitNpit(ind) = obj.ConConduitNpit(ind)+...
                            obj.ICConnections(i).Npit;
                        obj.ConConduitNpore(ind) = obj.ConConduitNpore(ind)+...
                            obj.ICConnections(i).Npore;
                    end
                end
            end
            obj.ConConduits = obj.ConConduits(1:count-1);
            obj.ConConduitASP = obj.ConConduitASP(1:count-1);
            obj.ConConduitNpit = obj.ConConduitNpit(1:count-1);
            obj.ConConduitNpore = obj.ConConduitNpore(1:count-1);
        end
        
    end
    
    
    methods (Static)
        function z = zeros(varargin)
            if (nargin == 0)
                % For zeros('Conduit')
                z = Conduit;
            elseif any([varargin{:}] <= 0)
                % For zeros with any dimension <= 0
                z = Conduit.empty(varargin{:});
            else
                % For zeros(m,n,...,'Conduit')
                % Use property default values
                z = repmat(Conduit,varargin{:});
            end
        end
    end
end