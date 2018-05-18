classdef Cluster < handle
    
    properties
        Conduits
        Rows
        ICConnections %Inter-conduit membranes
        Columns
        Tag %useful to speed up clustering
    end
    
    properties (Dependent)
        Connected
        NCCoduits %Non-Conducting Conduits
        NCConnections %Non-Conducting ICCs
        EConduits %Embolized Conduits
        NonConducting %Logical that makes sure embolisms don't completely block this cluster
    end
    
    methods
        
        function obj = Cluster(Tag,Rows,conduits)
            if nargin==0
                obj.Rows=0;
                obj.Conduits = Conduit.empty();
                obj.ICConnections = ICC.empty();
                obj.Tag = 0;
            end
            if nargin>0
                obj.Tag = Tag;
                obj.Rows=Rows;
                obj.Conduits=conduits;
%                 obj.Conduits(1).ClusterPointer = ClusterPointer(obj);
%                 for i=2:length(conduits)
%                     obj.Conduits(i).ClusterPointer = obj.Conduits(1).ClusterPointer;
%                 end
                for i=1:length(conduits)
                    conduits(i).Cluster=obj;
                end
                obj.ICConnections = ICC.empty();
            end
        end
        
        function AddConduits(obj,conduits)
            obj.Conduits=[obj.Conduits conduits];
            for i=1:length(conduits)
                conduits(i).Cluster=obj;
            end
%             conduits(1).ClusterPointer.Cluster = obj;
        end
        
        function IntegrateInto(obj,cluster)
            cluster.Conduits=[cluster.Conduits obj.Conduits];
%             cluster.ICConnections=[cluster.ICConnections obj.ICConnections];
            for i=1:length(obj.Conduits)
                obj.Conduits(i).Cluster=cluster;
            end

            delete(obj);
        end
        
        function Connected = get.Connected(obj)
            minimum=obj.Rows;
            maximum=0;
            for i=1:length(obj.Conduits)
                if max(obj.Conduits(i).Nodes(:,1))>maximum
                    maximum=max(obj.Conduits(i).Nodes(:,1));
                end
                if min(obj.Conduits(i).Nodes(:,1))<minimum
                    minimum=min(obj.Conduits(i).Nodes(:,1));
                end
            end
            if maximum==obj.Rows&&minimum==1
                Connected=true;
            else
                Connected=false;
            end
        end
        
        function updateColumns(obj)
            minimum=10000;
            maximum=0;
            for i=1:length(obj.Conduits)
                if obj.Conduits(i).Nodes(1,2)>maximum
                    maximum=obj.Conduits(i).Nodes(1,2);
                end
                if obj.Conduits(i).Nodes(1,2)<minimum
                    minimum=obj.Conduits(i).Nodes(1,2);
                end
            end
            obj.Columns=minimum:maximum;
        end
        
        function Isolate(obj,index)
            %This method isolates conduits that have no further connections
            %to a conductive cluster due to an embolism elsewhere
            if ~obj.Conduits(index).Embolized
               error('Conduit not embolized'); 
            end
            
            NEmb=~[obj.Conduits(index).ConConduits.Embolized];
            ConCon=obj.Conduits(index).ConConduits(NEmb);
            for i=1:length(ConCon)
                starts=false;
                ends=false; %Logicals that indicate whether section is isolated
                %If already passed over this connected conduit then skip it
                %and do the following
                if ConCon(i).NonConducting
                    continue
                end
                %Add embolized and connected vessels to visited  
                visited = Conduit.zeros(1,length(obj.Conduits));
                visitedInd = zeros(1,length(obj.Conduits));
                visited(1:2) = [obj.Conduits(index) ConCon(i)];
                %Index of first CE of every conduit for faster code when we
                %call ismember later
                visitedInd(1:2) = [obj.Conduits(index).OINodes(1),...
                    ConCon(i).OINodes(1)];
                visitCount = 2;
                if ConCon(i).isInlet
                    starts=true;
                end
                if ConCon(i).isOutlet
                    ends=true;
                end
                %Check if visiting an embolized conduit and disallow if so
                NEmb=~[ConCon(i).ConConduits.Embolized];
                %Also check whether conduit has already been visited
                visits = ConCon(i).ConConduits(...
                    and(~ismember(ConCon(i).ConConduits,...
                    visited(1:visitCount)),...
                    NEmb));
                %Keep on doing visits until no other connected conduit is
                %available
                while ~isempty(visits)

                    visited(visitCount+1:visitCount+length(visits)) =...
                        visits;
                    visitedInd(visitCount+1:visitCount+length(visits)) =...
                        [visits.firstOINode];
                    visitCount = visitCount + length(visits);

                    if  ~starts && any([visits.isInlet])                       
                        starts=true;
                    end
                    if  ~ends && any([visits.isOutlet])                    
                        ends=true;
                    end
                    CC=unique([visits.ConConduits]);
                    NEmb=~[CC.Embolized];%Not embolized logical
                    
                    %Check if connected conduit has been visited
                    %already from a different path                    
                    visits = CC(and(~ismember([CC.firstOINode],...
                        visitedInd(1:visitCount)),...
                        NEmb));
                end
                %Isolated cluster or not?
                if ~starts || ~ends
                   arrayfun(@(x) x.Isolate,visited(1:visitCount))
                end
            end 
        end
        
        function NonConducting = get.NonConducting(obj)
            found=false;
            for i=1:length(obj.Conduits)
                if ~obj.Conduits(i).NonConducting
                    found=true;
                    break
                end
            end
            if found
                NonConducting=false;
            else
                NonConducting=true;
            end
        end      
    end
    
    methods (Static)
        function z = zeros(varargin)
            if (nargin == 0)
                % For zeros('Cluster')
                z = Cluster;
            elseif any([varargin{:}] <= 0)
                % For zeros with any dimension <= 0
                z = Cluster.empty(varargin{:});
            else
                % For zeros(m,n,...,'Cluster')
                % Use property default values
                z = repmat(Cluster,varargin{:});
            end
        end
    end
end