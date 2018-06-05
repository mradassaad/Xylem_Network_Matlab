classdef XylemNet < handle
    
    properties
        Conduits %Functional Conduits
        ConSorted
        ICConnections %Inter-Conduit Connections
        Clusters
        Kend %Vector containing last row of conductances
        %The following properties are useful for stitching
        Aend %End section vessel areas
        CBProw
        CBPcol
        pickedConx
    end
    
    properties (SetAccess = immutable)
        Pc
        Pe %Probability of node becoming vessel end
        NPe %Probability of node becoming vessel beginning
        Size %dimension of the grid on which Conduits were created
        ConduitScheme
        EndWallScheme
        Dcross %diameters at 4 cross sections useful to compute resistivity
    end
    
    properties (Dependent)
        Embolized %Check if none of the connected clusters are conducting anymore
    end
    
    properties (Constant)
        gamma = 0.072; %Pa.m
    end
    
    methods
        
        function obj = XylemNet(row,column,Pe,NPe,Pc,Lce,Dc,Dc_cv,varargin)
            %Class constructor
            
            %Case where there is no input
            if nargin == 0
                obj.Conduits = Conduit.empty;
                obj.Size = 0;
                obj.ICConnections = EndWall.empty;
                
            elseif nargin > 0
                %Parse required inputs
                
                %expectedConduitSchemes:
                %1 - Random diameter assignment
                %2 - Assign diameters in ranks
                %expectedEndWallSchemes:
                %1 - Takes as inputs Dp,Dp_cv,Fc,Fpf
                reqInputs = inputParser;
                validPosNum = @(x) isnumeric(x) && all(x > 0);
                addRequired(reqInputs,'row',@(x) validPosNum(x));
                addRequired(reqInputs,'column',@(x) validPosNum(x));
                addRequired(reqInputs,'Pe',@(x) validPosNum(x) ...
                    && x<=1);
                addRequired(reqInputs,'NPe',@(x) validPosNum(x) ...
                    && x<=1);
                addRequired(reqInputs,'Pc',@(x) validPosNum(x) ...
                    && isscalar(x) && x<=1);
                addRequired(reqInputs,'Lce',@(x) validPosNum(x)); %meters
                addRequired(reqInputs,'Dc',@(x) validPosNum(x) ...
                    && isscalar(x) && x<200e-6 && x>1e-6); %meters
                addRequired(reqInputs,'Dc_cv',@(x) validPosNum(x) ...
                    && isscalar(x));
                parse(reqInputs,row,column,Pe,NPe,Pc,Lce,Dc,Dc_cv);
                
                obj.Size = [row column];
                obj.Pc = Pc;
                obj.Pe = Pe;
                obj.NPe = NPe;
                
                %Parse optional inputs
                optInputs = inputParser;
                addRequired(optInputs,'Dp',@(x) validPosNum(x) ...
                    && isscalar(x) && x<200e-9 && x>1e-10); %meters
                addRequired(optInputs,'Dm',@(x) validPosNum(x) ...
                    && isscalar(x) && x<200e-6 && x>1e-7); %meters
                addRequired(optInputs,'k_ASP',@(x) validPosNum(x) ...
                    && isscalar(x));
                addRequired(optInputs,'lam_ASP',@(x) validPosNum(x) ...
                    && isscalar(x));
                addRequired(optInputs,'Fc',@(x) validPosNum(x) ...
                    && isscalar(x) && x<=1);
                addRequired(optInputs,'Fpf',@(x) validPosNum(x) ...
                    && isscalar(x) && x<=1);
                addRequired(optInputs,'Fap',@(x) validPosNum(x) ...
                    && isscalar(x) && x<=1);
                addRequired(optInputs,'e_mean',@(x) validPosNum(x) ...
                    && isscalar(x));
                addRequired(optInputs,'e_cv',@(x) validPosNum(x) ...
                    && isscalar(x));
                addRequired(optInputs,'Tm',@(x) validPosNum(x) ...
                    && isscalar(x)); %meters
                addRequired(optInputs,'Lp',@(x) validPosNum(x) ...
                    && isscalar(x)); %meters
                addRequired(optInputs,'ASPcalcmethod',@(x) ischar(x));
                parse(optInputs,varargin{:});
                
                Dp = optInputs.Results.Dp;
                Dm = optInputs.Results.Dm;
                k_ASP = optInputs.Results.k_ASP;
                lam_ASP = optInputs.Results.lam_ASP;
                Fc = optInputs.Results.Fc;
                Fpf = optInputs.Results.Fpf;
                Fap = optInputs.Results.Fap;
                e_mean = optInputs.Results.e_mean;
                e_cv = optInputs.Results.e_cv;
                Tm = optInputs.Results.Tm;
                Lp = optInputs.Results.Lp;
                ASPcalcmethod = optInputs.Results.ASPcalcmethod;
                %Generate Conduits with constant conduit element length
                [obj.Conduits,obj.ConSorted,obj.CBProw,obj.CBPcol] =...
                    generateConduits(obj.Size,Lce,obj.Pe,obj.NPe);
                
                %Generate ICCs and add to corresponding conduits
                %First determine potential connections defined as any two
                %adjacent nodes that are part of different conduits 
                potConx = findConx(obj.Conduits,obj.Size(1));
                %Choose potential connections randomly depending on
                %prescribed probability Pc
                pickedConx = pickConx(potConx,obj.Pc);
                %Save that array in case we need to create a new xylem
                %network with the same ICC positions
                obj.pickedConx = pickedConx;
                %With the chosen adjacent nodes, create ICC objects
                [obj.ICConnections,obj.Clusters] = ...
                    generateICCs(obj,Dp,Dm,k_ASP,lam_ASP,Fc,Fpf,Fap,Tm,Lp,...
                    ASPcalcmethod);
                %updateConduits is a XylemNet method that deletes all
                %conduits that aren't part of any cluster and are therefore
                %not useful for the overall hydraulic pathway of the
                %segment
                obj.updateConduits();
                if isempty(obj.ICConnections)
                    disp('No Network')
                    return
                end
                
                Dstd=Dc_cv*Dc; %m
                Dm=log(Dc^2/sqrt(Dstd^2+Dc^2));
                Ds=sqrt(log(Dstd^2/(Dc^2)+1));
                Dcs=lognrnd(Dm,Ds,1,length(obj.Conduits));
                Dcs=sort(Dcs);
                tempLen=[obj.Conduits.Length];
                [~,ILen]=sort(tempLen);
                fn = DcMapping(length(obj.Conduits),...
                    floor(length(obj.Conduits)/20));
                for k=1:length(obj.Conduits)
                    obj.Conduits(k).updateDiameter(...
                        Dcs(fn(ILen==k)));
                end
                
                for i=1:length(obj.ConSorted)
                    if isa(obj.ConSorted{i},'Conduit')
                        obj.ConSorted{i}=obj.ConSorted{i}(isvalid(obj.ConSorted{i}));
                    end
                end
                for i=1:length(obj.Conduits)
                    obj.Conduits(i).addConConduits;
                    obj.Conduits(i).updateConduit;
                end

                for i=1:length(obj.ICConnections)
                    obj.ICConnections(i).updateMeanArea;
                end
                
                es = sort(abs(normrnd(e_mean,e_cv*e_mean,...
                    1,length(obj.ICConnections))));
                [~,I] = sort([obj.ICConnections.Am]);
                for i=1:length(obj.ICConnections)
                    obj.ICConnections(i).computeKmASP(es(I==i));
                end
                
                for i=1:length(obj.Clusters)
                    obj.Clusters(i).updateColumns;
                end
                for i=1:length(obj.Conduits)
                    obj.Conduits(i).addConConduitASP;
                end
                
                x=obj.Size(1):obj.Size(1):(obj.Size(2))*obj.Size(1);
                obj.Kend=zeros(1,obj.Size(2));
                obj.Aend=zeros(1,obj.Size(2));
                for i=1:length(obj.Kend)
                    temp=obj.getConduits(x(i));
                    if ~isempty(temp)
                        obj.Kend(i)=temp.Kce;
                        obj.Aend(i) = pi.*(temp.Diameter./2).^2;
                    end
                end
                firstOINodes = [obj.Conduits.firstOINode];
                lastOINodes = [obj.Conduits.lastOINode];
                OINodes = [firstOINodes;lastOINodes];
                OINodes = OINodes(:);
                [obj.CBProw,obj.CBPcol]=ind2sub(obj.Size,OINodes);
                if length(obj.ICConnections)>1
                    obj.pickedConx={obj.ICConnections.Nodes};
                end
                
                %Cross sectional area of conduits at 4 cross sections
                obj.Dcross =...
                    [obj.getD(1); obj.getD(2); obj.getD(3); obj.getD(4)];
            end
            
        end
        
        %Network properties with embolization
        function [K,Fo,Embolized,Penetrated]=AirSeed(obj,Pair,ClusNb,EmbNb,SorL,RankSize,VorLPP)
            %EmbNb being the number of random embolizations, in the
            %beginning, for every cluster
            %ClusNb determines the clusters to be embolized
            getPflag=false;
            if EmbNb>0
                getPflag=true;
                obj.RepairAll;
                for k=ClusNb
                    if ~exist('SorL','var')
                        available=obj.Clusters(k).Conduits;
                    elseif strcmp(VorLPP,'V')
                        Dias = [obj.Clusters(k).Conduits.Diameter];
                        [~,I] = sort(Dias);
                        if strcmp(SorL,'S')
                            available=obj.Clusters(k).Conduits(I(1:RankSize));
                        elseif strcmp(SorL,'M')
                            available=obj.Clusters(k).Conduits(I(floor((end-RankSize)/2):ceil((end+RankSize)/2)));
                        elseif strcmp(SorL,'L')
                            available=obj.Clusters(k).Conduits(I(end-RankSize+1:end));
                        end
                    elseif strcmp(VorLPP,'LPP')
                        Dp = [obj.Clusters(k).ICConnections.Dp];
                        [~,I] = sort(Dp);
                        count=0;
                        i=1;
                        available = Conduit.empty;
                        while count < RankSize
                            if strcmp(SorL,'S')
                                Con1 = obj.Clusters(k).ICConnections(I(i)).ConConduits(1);
                                Con2 = obj.Clusters(k).ICConnections(I(i)).ConConduits(2);
                            elseif strcmp(SorL,'L')
                                Con1 = obj.Clusters(k).ICConnections(I(end-i+1)).ConConduits(1);
                                Con2 = obj.Clusters(k).ICConnections(I(end-i+1)).ConConduits(2);
                            end
                            if ~ismember(Con1,available)
                                available(count+1) = Con1;
                                count = count+1;
                            end
                            if ~ismember(Con2,available)
                                available(count+1) = Con2;
                                count = count+1;
                            end
                            i = i+1;
                        end
                    end
                    for i=1:EmbNb
                        if ~isempty(available)
                            pick=round(rand*(length(available)-1))+1;
                            obj.Embolize(available(pick));
                            available(pick)=[];
                        end
                    end
                end
            end
            %             Embolized = cell(1,100);
            %             Penetrated = cell(1,100);
            [CProp,Penetrated]=obj.PropagateAirSeed(Pair);
            Embolized = CProp;
            %             count = 2;
            
            while ~isempty(CProp)
                getPflag=true;
                for i=1:length(CProp)
                    obj.Embolize(CProp(i))
                end
                [CProp,Pen] = obj.PropagateAirSeed(Pair);
                Embolized = [Embolized CProp]; %#ok<AGROW>
                Penetrated = [Penetrated Pen]; %#ok<AGROW>
                %                 count = count + 1;
                
            end
            Embolized = unique(Embolized);
            Penetrated = unique(Penetrated);
            if getPflag
                [~,Fo,K,~,~] = compute_hydraulics(obj,0.1,0);
            else
                K=-1;Fo=-1;
            end
            %             Embolized = Embolized(1:count-2);
            %             Penetrated = Penetrated(1:count-2);
        end
        
        
        
        function [GCE, GCEXdata, GCEYdata, GConHyd, GConCav, GConIn...
                , GConOut] = createGraphs(obj)
            temp = [obj.Conduits];
            ACEsize = sum([temp.CECount]+1);
            ACE = zeros(ACEsize);
            lastpos = 1;
            Xdata = zeros(1,ACEsize);
            Ydata = zeros(1,ACEsize);
            %Conduit Graphs
            AConsize = length(temp);
            AConHyd = zeros(AConsize);
            AConCav = zeros(AConsize);
            in = 0;
            incount = 1;
            out = 0;
            outcount = 1;
            lastposCon = 1;
            for k=1:length(obj.Clusters)
                if obj.Clusters(k).NonConducting
                    continue
                end
                nrows=obj.Size(1);
                ncolumns=length(obj.Clusters(k).Columns);
                firstColumn=obj.Clusters(k).Columns(1);
                un=false(1,nrows*ncolumns);
                
                %Fill in values where n~=m and are part of a conduit
                
                cc=obj.Clusters(k).Conduits(~[obj.Clusters(k).Conduits.NonConducting]);
                [~,I]=sort([cc.firstOINode]);
                cc=cc(I);
                AAHyd = zeros(length(cc));
                AACav = zeros(length(cc));
                NN={cc.OINodes};
                lastpos2 = lastpos;
                objsize = [obj.Size(1) obj.Size(2)];
                for i = 1:length(NN)
                    [x,y] = ind2sub(objsize,NN{i});
                    Xdata(lastpos2:length(NN{i})+lastpos2-1) = y;
                    Ydata(lastpos2:length(NN{i})+lastpos2-1) = x;
                    lastpos2 = length(NN{i})+lastpos2;
                end
                AA=zeros(sum([cc.CECount]+1));
                pos=[0 cumsum( [cc.CECount]+1)];
                for i=1:length(NN)
                    UNodes = NN{i}-nrows.*(firstColumn-1);
                    un(UNodes)=true;
                    OINodes=pos(i)+1:pos(i+1);
                    
                    AA(size(AA,1)*(OINodes(2:end)-1)+OINodes(1:end-1))=cc(i).Kce; %mm^3/Pa day
                    AA(size(AA,1)*(OINodes(1:end-1)-1)+OINodes(2:end))=cc(i).Kce; %mm^3/Pa day
                end
                unun=find(un);
                
                ConCon=obj.Clusters(k).ICConnections(~[obj.Clusters(k).ICConnections.NonConducting]);
                NodeNode={ConCon.OINodes};
                
                for i=1:length(NodeNode)
                    pitOINodes = NodeNode{i}-nrows.*(firstColumn-1);
                    pitOINodes(1) = find(unun==pitOINodes(1));
                    pitOINodes(2) = find(unun==pitOINodes(2));
                    
                    AA(pitOINodes(2),pitOINodes(1))=ConCon(i).Km; %mm^3/Pa day
                    AA(pitOINodes(1),pitOINodes(2))=ConCon(i).Km; %mm^3/Pa day
                end
                
                %Conduit Graphs
                for i=1:length(cc)
                    temp = cc(i).ICConnections;
                    for j = 1:length(temp)
                        if isequal(temp(j).ConConduits(1),cc(i))
                            [~,ind] = ismember(temp(j).ConConduits(2),cc(i+1:end));
                            AAHyd(i,ind+i) = AAHyd(i,ind+i)+temp(j).Km;
                            AAHyd(ind+i,i) = AAHyd(i,ind+i);
                            AACav(i,ind+i) = max(AACav(i,ind+i),(4*0.072*10^(-6))/temp(j).Dp); %MPa
                            AACav(ind+i,i) = AACav(i,ind+i);
                        end
                    end
                    if cc(i).isinput
                        in(incount) = i+lastposCon-1; %#ok<AGROW>
                        incount = incount + 1;
                    elseif cc(i).isoutput(obj.Size)
                        out(outcount) = i+lastposCon-1; %#ok<AGROW>
                        outcount = outcount + 1;
                    end
                end
                
                ACE(lastpos:sum([cc.CECount]+1)+(lastpos-1),lastpos:sum([cc.CECount]+1)+(lastpos-1))=AA;
                lastpos = sum([cc.CECount]+1)+lastpos;
                %Conduit Graphs
                AConHyd(lastposCon:length(cc)+(lastposCon-1),lastposCon:length(cc)+(lastposCon-1))=AAHyd;
                AConCav(lastposCon:length(cc)+(lastposCon-1),lastposCon:length(cc)+(lastposCon-1))=AACav;
                lastposCon = length(cc)+lastposCon;
            end
            GCE = graph(ACE);
            GCEXdata = Xdata;
            GCEYdata = Ydata;
            GConHyd = graph(AConHyd);
            GConCav = graph(AConCav);
            GConIn = in;
            GConOut = out;
        end
        
        function Conduits = updateConduits(obj)
            %This function removes non functional conduits, non connecting
            %clusted, and non functional ICCs
            Conduits = obj.Conduits([obj.Conduits.Functional]);
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
%             obj.plotNet    
            prune(obj)  
        end
        
        function [Pmatrix,singular] = getPressures(obj,Pin,Pout)
            %This function returns the matrix of pressures and the matrix
            %of conductances
            if obj.Embolized
                %                 warning('Completely embolized')
                Pmatrix = zeros(obj.Size);
                singular=0;
                return
            end
            %singular is a logical that checks whether all isolated
            %conduits are removed
            singular=false;
            %Initialize pressure matrix
            Pmatrix = zeros(obj.Size);
            
            %Solve by cluster: works since clusters are independent
            %hydraulically
            for k=1:length(obj.Clusters)
                %Check if cluster connects segment inlet to outlet
                if obj.Clusters(k).NonConducting
                    continue
                end
                nrows=obj.Size(1);
                ncolumns=length(obj.Clusters(k).Columns);
                firstColumn=obj.Clusters(k).Columns(1);
                %un is a bolean array. It is true where the linear index of
                %the segment grid contains a conduit element or ICC. It is
                %used to reduce the sizes of A and B to reduce memory use
                %and computational cost. Later, un will be used to compute
                %unun which contains the linear indeces of conduit elements
                %and ICCs in A after A has been reduced to include only
                %those nodes that aren't empty.
                un=false(1,nrows*ncolumns);
                %PP is the pressure matrix of this cluster
                PP = zeros([nrows ncolumns]);
                
                %B contains the pressure at the inlet and outlets
                B=zeros(nrows*ncolumns,1);
                
                %Identify inlet and outlet nodes in the new numbering
                %system
                inlets = 1:nrows:nrows*(ncolumns-1)+1;
                outlets = nrows:nrows:nrows*(ncolumns-1)+nrows;
                inOut = [inlets, outlets]';
                B(inOut(1:length(inOut)/2))=Pin;
                B(inOut(length(inOut)/2+1:end))=Pout;
                
                %Fill in values where n~=m and are part of a conduit
                %cc contains all conducting conduits of this cluster (not
                %embolized nor isolated)
                cc = obj.Clusters(k).Conduits(...
                    ~[obj.Clusters(k).Conduits.NonConducting]);
                %Make sure we study conduits consecutively depending on
                %their position in the segment
                [~,I]=sort([cc.firstOINode]);
                cc=cc(I);
                %NN contains the one index position of the consecutive
                %conduits. One index positions are the linear indeces.
                NN={cc.OINodes};
                
                %A will contain the conductances
                A=zeros(sum([cc.CECount]+1));
                pos=[0 cumsum( [cc.CECount]+1)];
                for i=1:length(NN)
                    %Correction of linear indeces in case the first cluster
                    %column is not column 1 of the segment. This is because
                    %matrix A needs to start from the first column
                    %otherwise a lot of memory will be used up.
                    UNodes = NN{i}-nrows.*(firstColumn-1);
                    un(UNodes)=true;
                    OINodes=pos(i)+1:pos(i+1);
                    %Identify whether conduit corresponding to these nodes
                    %is an inlet. If it is then A(OINodes(1),OINodes(1))
                    %should equal 1.
                    if mod(UNodes(1)-1,obj.Size(1))==0
                        A(OINodes(2),OINodes(1))=-cc(i).Kce; %m^3/MPa.s
                        %Only reason this is -1 is because the sum step
                        %below reverses it to 1
                        A(OINodes(1),OINodes(1))=-1;
                        OINodes=OINodes(2:end);
                    %Identify whether conduit corresponding to these nodes
                    %is an outlet. If it is then A(OINodes(1),OINodes(1))
                    %should equal 1.
                    elseif mod(UNodes(end),obj.Size(1))==0
                        A(OINodes(end-1),OINodes(end))=-cc(i).Kce; %m^3/MPa.s
                        A(OINodes(end),OINodes(end))=-1;
                        OINodes=OINodes(1:end-1);
                    end
                    %For the rest of the nodes of inlet and outlet
                    %conduits, or for all the nodes of all other conduits,
                    %insert -Kce on both sides of the main diagonal of A
                    A(size(A,1)*(OINodes(2:end)-1)+OINodes(1:end-1)) =...
                        -cc(i).Kce; %m^3/MPa.s
                    A(size(A,1)*(OINodes(1:end-1)-1)+OINodes(2:end)) =...
                        -cc(i).Kce; %m^3/MPa.s
                end
                unun=find(un);
                
                ICCs = obj.Clusters(k).ICConnections(...
                    ~[obj.Clusters(k).ICConnections.NonConducting]);
                ICCind={ICCs.OINodes};
                
                %Add ICC conductances to matrix A
                for i=1:length(ICCind)
                    %Same as before correct the linear index nodes of ICCs
                    pitOINodes = ICCind{i}-nrows.*(firstColumn-1);
                    %Find the position in A where -Km should be. Because A
                    %contains only nodes that have a conduit element or an
                    %ICC in them, its indeces do not correspond to those of
                    %the segment grid. Therefore, find(unun==pitOINodes(1))
                    %finds the indeces of the ICC in the reduced A.
                    pitOINodes(1) = find(unun==pitOINodes(1));
                    pitOINodes(2) = find(unun==pitOINodes(2));
                    
                    A(pitOINodes(2),pitOINodes(1))=-ICCs(i).Km; %m^3/MPa.s
                    A(pitOINodes(1),pitOINodes(2))=-ICCs(i).Km; %m^3/MPa.s
                end
                %Now B is also reduced like A to contain only nodes that
                %contain a conduit element or ICC
                B=B(un);
                
                %Fill in values where n==m, The diagonal should contain the
                %sum of all conductances going through a certain node
                %except for the inlets and outlets.
                A(eye(size(A),'like',true))= -sum(A,2);
                
                A=sparse(A);
                B=sparse(B);
                
                try
                    lastwarn('');
                    PP(un)= A\B;
                    error(lastwarn)
                catch
                    singular=true;
                    disp('Singular')
                end
                
                Pmatrix(:,obj.Clusters(k).Columns) =...
                    Pmatrix(:,obj.Clusters(k).Columns) + PP;
            end
            
        end
        
        function Embolize(obj,input)
            if isa(input,'Conduit')
                input.Embolize;
                for j=1:length(input.Cluster.Conduits)
                    if isequal(input.Cluster.Conduits(j),input)
                        break
                    end
                end
                input.Cluster.Isolate(j);
                
            elseif input>length(obj.Conduits)
                error('input greater than total number of functional conduits')
            else
                found=false;
                for i=1:length(obj.Clusters)
                    for j=1:length(obj.Clusters(i).Conduits)
                        if isequal(obj.Clusters(i).Conduits(j),obj.Conduits(input))
                            found=true;
                            break
                        end
                    end
                    if found
                        break
                    end
                end
                obj.Conduits(input).Embolize;
                obj.Clusters(i).Isolate(j);
            end
            
        end
        
        function [Fv,Fh]= getFlow(obj,Pmatrix)
            %This function returns the flowfield in the whole grid.
            %There will be two outputs to this function: the vertical and
            %horizontal flow
            %The convention will be to take up and right as positive (down
            %and left as negative)
            Fv=zeros(obj.Size(1),obj.Size(2));
            Fh=zeros(obj.Size(1),obj.Size(2));
            for k=1:length(obj.Clusters)
                if obj.Clusters(k).NonConducting
                    continue
                end
                nrows=obj.Size(1);
                ncolumns=length(obj.Clusters(k).Columns);
                firstColumn=obj.Clusters(k).Columns(1);
                
                Fvv=zeros(nrows,ncolumns);
                Fhh=zeros(nrows,ncolumns);
                
                PP=Pmatrix(:,obj.Clusters(k).Columns);
                %Fill in values for Kv (To clarify, the values in this matrix
                %will only include flows from a node to an upper node. F34 will
                %give the flow from node 3 to 4 if it is positive. From 4 to 3
                %if it is negative)
                cc=obj.Clusters(k).Conduits(~[obj.Clusters(k).Conduits.NonConducting]);
                NN={cc.OINodes};
                for i=1:length(NN)
                    OINodes = NN{i}-nrows.*(firstColumn-1);
                    OINodes=OINodes(1:end-1);
                    
                    [I,J] = ind2sub(size(Fvv),OINodes);
                    Fvv(I,J) = cc(i).Kce*...
                        (PP(I,J) - PP(I+1,J));
                end
                
                %Fill in values for Kh (To clarify, the values in this matrix
                %will only include flows from a node to a another node on the right.
                %F39 will give the flow from node 3 to 9 if it is positive. From 9 to 3
                %if it is negative)
                
                ConCon=obj.Clusters(k).ICConnections(~[obj.Clusters(k).ICConnections.NonConducting]);
                NodeNode={ConCon.OINodes};
                
                for i=1:length(NodeNode)
                    pitOINodes = NodeNode{i}-nrows.*(firstColumn-1);
                    if pitOINodes(1)<pitOINodes(2)
                        [I,J] = ind2sub(size(Fhh),pitOINodes);
                        Fhh(I(1),J(1)) = ConCon(i).Km*...
                            (PP(I(1),J(1)) - PP(I(2),J(2)));
                    else
                        [I,J] = ind2sub(size(Fhh),pitOINodes);
                        Fhh(I(2),J(2)) = ConCon(i).Km*...
                            (PP(I(2),J(2)) - PP(I(1),J(1)));
                    end
                end
                
                Fv(:,obj.Clusters(k).Columns) = ...
                    Fv(:,obj.Clusters(k).Columns)+Fvv; % m3/s
                Fh(:,obj.Clusters(k).Columns) =...
                    Fh(:,obj.Clusters(k).Columns)+Fhh;
            end
        end
        
        function [Infected,EW_penetrated]=PropagateAirSeedbackup(obj,Pair)
            
            embConduits = obj.Conduits([obj.Conduits.Embolized]);
            Infected = Conduit.zeros(1,1000);
            InfectedN1 = zeros(1,1000); %First index of infected conduits
            EW_penetrated = EndWall.empty;
            count=1;
            Pcomp=0.1; %atmospheric pressure (MPa)
            for i=1:length(embConduits)
                for j=1:length(embConduits(i).ICConnections)
                    check=0;
                    if ~embConduits(i).ICConnections(j).ConConduits(1).NonConducting...
                            && ~ismember(embConduits(i).ICConnections(j).ConConduits(1).firstOINode,InfectedN1(1:count-1))
                        check=1;
                    elseif ~embConduits(i).ICConnections(j).ConConduits(2).NonConducting...
                            && ~ismember(embConduits(i).ICConnections(j).ConConduits(2).firstOINode,InfectedN1(1:count-1))
                        check=2;
                    end
                    %difference between bubble pressure and adjacent node
                    difference=Pair-embConduits(i).ICConnections(j).ASP;
                    if check==1
                        if Pcomp <= difference
                            Infected(count)=embConduits(i).ICConnections(j).ConConduits(1);
                            InfectedN1(count)=embConduits(i).ICConnections(j).ConConduits(1).firstOINode;
                            EW_penetrated(count) = embConduits(i).ICConnections(j);
                            count=count+1;
                        end
                    elseif check==2
                        if Pcomp <= difference
                            Infected(count)=embConduits(i).ICConnections(j).ConConduits(2);
                            InfectedN1(count)=embConduits(i).ICConnections(j).ConConduits(2).firstOINode;
                            EW_penetrated(count) = embConduits(i).ICConnections(j);
                            count=count+1;
                        end
                    end
                end
            end
            Infected = Infected(1:count-1);
            EW_penetrated = EW_penetrated(1:count-1);
        end
        
        function [Infected,EW_penetrated]=PropagateAirSeed(obj,Pair)
            
            embConduits = obj.Conduits([obj.Conduits.Embolized]);
            Infected = Conduit.zeros(1,1000);
            InfectedN1 = zeros(1,1000); %First index of infected conduits
            EW_penetrated = zeros(1,1000);
            count=1;
            Pcomp=0.1; %atmospheric pressure (MPa)
            for i=1:length(embConduits)
                ConductiveBool =...
                    ~[embConduits(i).ConConduits.NonConducting];
                Conductive = embConduits(i).ConConduits(ConductiveBool);
                ConductiveASP = embConduits(i).ConConduitASP(ConductiveBool);
                for j=1:length(Conductive)
                    difference=Pair-ConductiveASP(j);
                    if Pcomp <= difference
                        Infected(count) = Conductive(j);
                        embConduits(i).ConConduits(j).firstOINode;
                        InfectedN1(count) = Conductive(j).firstOINode;
                        embConduits(i).ConConduitASP(j);
                        EW_penetrated(count) =...
                            ConductiveASP(j);
                        count=count+1;
                    end

                end
            end
            Infected = Infected(1:count-1);
            EW_penetrated = EW_penetrated(1:count-1);
        end
        
        %Auxiliary functions
        function [CD,CL,LPP,NPM,NCC,TNC]=getProperties(obj)
            %               ? Conduit Diameter
            %               ? Conduit Length
            %             	? Largest Pit Pore
            %             	? Number of PitMemranes
            %             	? Number of connected conduits
            %             		(1) Mean
            %             		(2) Std
            %             		(3) Histogram
            %             	? Total Number of conduits
            CD(1)=mean([obj.Conduits.Diameter]); %m
            CD(2)=std([obj.Conduits.Diameter]); %m
            
            CL(1)=mean([obj.Conduits.Length]); %m
            CL(2)=std([obj.Conduits.Length]); %m
            
            LPP(1)=mean([obj.ICConnections.Dp]);
            LPP(2)=std([obj.ICConnections.Dp]);
            
            NPM(1)=mean([obj.Conduits.ICCCount]);
            NPM(2)=std([obj.Conduits.ICCCount]);
            
            NCC(1)=mean([obj.Conduits.NConConduits]);
            NCC(2)=std([obj.Conduits.NConConduits]);
            
            TNC=length(obj.Conduits);
        end
        
%         function plot(obj)
%             if isempty(obj.GCE)
%                 disp('Initializing graph')
%                 [obj.GCE, obj.GCEXdata, obj.GCEYdata,...
%                     obj.GConHyd, obj.GConCav, obj.GConIn, obj.GConOut]...
%                     = createGraphs(obj);
%             end
%             figure;plot(obj.GCE,'Marker','none',...
%                 'XData',obj.GCEXdata,...
%                 'YData',obj.GCEYdata)
%         end
%         
%         function plotCon(obj)
%             if isempty(obj.GCE)
%                 disp('Initializing graph')
%                 [obj.GCE, obj.GCEXdata, obj.GCEYdata,...
%                     obj.GConHyd, obj.GConCav, obj.GConIn, obj.GConOut]...
%                     = createGraphs(obj);
%             end
%             figure;
%             pConHyd = plot(obj.GConHyd,'layout','layered',...
%                 'Sources',obj.GConIn,'Sinks',obj.GConOut);
%             highlight(pConHyd,obj.GConIn,'NodeColor','g','MarkerSize',5)
%             highlight(pConHyd,obj.GConOut,'NodeColor','r','MarkerSize',5)
%         end
        
        function visualize(obj)
            % create a default color map ranging from red to light pink for
            % LPP
            %             [~,edgesCD,~] = histcounts([obj.Conduits.Diameter].*10^6);
            %             lengthCD = length(edgesCD);
            %             widthCD = edgesCD(2)-edgesCD(1);
            %             shiftCD = floor(edgesCD(1)/widthCD);
            %             [~,edgesLPP,~] = histcounts([obj.ICConnections.Dpm].*10^9);
            %             lengthLPP = length(edgesLPP);
            %             widthLPP = edgesLPP(2)-edgesLPP(1);
            %             shiftLPP = floor(edgesLPP(1)/widthLPP);
            %Accounts for when the first edge is a multiple>2 larger than
            %the width
%             figure
            %             colors_LPP = colormap(cool(lengthLPP));
            %             colors_CD = colormap(gray(lengthCD));
            ax = axes;
            %             grid on
            for i=1:length(obj.Conduits)
                if ~obj.Conduits(i).Embolized
                    %                     lineColor = colors_CD(int8(lengthCD-floor((obj.Conduits(i).Diameter.*10^6)./widthCD)+shiftCD),:);
                    %                     plot(ax,obj.Conduits(i).Nodes(:,2),obj.Conduits(i).Nodes(:,1),'Color',lineColor,'LineWidth',1.3);
                    plot(ax,obj.Conduits(i).Nodes(:,2),obj.Conduits(i).Nodes(:,1),'k','LineWidth',2.5);
                    if obj.Conduits(i).NonConducting
                        plot(ax,obj.Conduits(i).Nodes(:,2),obj.Conduits(i).Nodes(:,1),'m','LineWidth',2.5);
                    end
                else
                    plot(ax,obj.Conduits(i).Nodes(:,2),obj.Conduits(i).Nodes(:,1),'c','LineWidth',2.5);
                end
                hold on
                %                 grid on
                
                for j=1:obj.Conduits(i).Nm
                    if ~obj.Conduits(i).ICConnections(j).NonConducting
                        %                         lineColor = colors_LPP(int8(lengthLPP-floor((obj.Conduits(i).ICConnections(j).Dpm.*10^9)./widthLPP)+shiftLPP),:);
                        %                         plot(ax,obj.Conduits(i).ICConnections(j).Nodes(:,2),obj.Conduits(i).ICConnections(j).Nodes(:,1),'Color',lineColor);
                        plot(ax,obj.Conduits(i).ICConnections(j).Nodes(:,2),obj.Conduits(i).ICConnections(j).Nodes(:,1),'r','LineWidth',1.3);
                    else
                        if obj.Conduits(i).ICConnections(j).Embolized
                            plot(ax,obj.Conduits(i).ICConnections(j).Nodes(:,2),obj.Conduits(i).ICConnections(j).Nodes(:,1),'c','LineWidth',1.3);
                        else
                            plot(ax,obj.Conduits(i).ICConnections(j).Nodes(:,2),obj.Conduits(i).ICConnections(j).Nodes(:,1),'m','LineWidth',1.3);
                        end
                    end
                    hold on
                end
            end
            ax.FontSize = 12;
            ax.LineWidth = 1.3;
            ax.FontWeight = 'bold';
            ax.XTick = [];
            ax.YTick = [];
            %             if ~isempty(obj.Conduits)
            %                 title(strcat('Vessel element length = ', num2str(obj.Conduits(1).CEs(1).Length.*10^3),' mm'))
            %             end
            xlim([0 obj.Size(2)+1]);
            ylim([1 obj.Size(1)]);
        end
        
        function cons = getConduits(obj,Nodes)
            %Returns Conduits corresponding to the input Nodes
            cons=Conduit.empty();
            if size(Nodes,2)==1
                tempNodes=[mod(Nodes,obj.Size(1)) floor(Nodes./obj.Size(1))+1];
                tempNodes(tempNodes(:,1)==0,2)=tempNodes(tempNodes(:,1)==0,2)-1;tempNodes(tempNodes(:,1)==0,1)=obj.Size(1);
                Nodes=tempNodes;
            end
            
            for j=1:size(Nodes,1)
                for i=1:length(obj.ConSorted{Nodes(j,2)})
                    NN=obj.ConSorted{Nodes(j,2)}(i).Nodes(:,1);
                    if Nodes(j,1)>=NN(1) && Nodes(j,1)<=NN(end)
                        cons(j) = obj.ConSorted{Nodes(j,2)}(i);
                    end
                end
            end
        end
        
        function [Rs] =	Resistivity(obj)
            %Total network conductivity
            Lce = obj.Conduits(1).CEs(1).Length; %m
            [~,~,Ktot,~] = compute_hydraulics(obj,0.1,0);
            Rtot = 1./Ktot; %MPa.s/m3
            Rc = Rtot./ (Lce*obj.Size(1)); %MPa.s/m4
            
            D1 = obj.Dcross(1,:);
            D2 = obj.Dcross(2,:);
            D3 = obj.Dcross(3,:);
            D4 = obj.Dcross(4,:);
            
            RC1 = Rc * length(D1(D1>0)); %MPa.s/m4
            RC2 = Rc * length(D2(D2>0)); %MPa.s/m4
            RC3 = Rc * length(D3(D3>0)); %MPa.s/m4
            RC4 = Rc * length(D4(D4>0)); %MPa.s/m4
            RC = mean([RC1 RC2 RC3 RC4]); %MPa.s/m4 average resistivity of one conduit
            
            RL1 = (sum(pi.*(D1(D1>0).^4)./(128 * 1.002e-3)))^(-1)*length(D1(D1>0))./10^6; %MPa.s/m4
            RL2 = (sum(pi.*(D2(D2>0).^4)./(128 * 1.002e-3)))^(-1)*length(D2(D2>0))./10^6; %MPa.s/m4
            RL3 = (sum(pi.*(D3(D3>0).^4)./(128 * 1.002e-3)))^(-1)*length(D3(D3>0))./10^6; %MPa.s/m4
            RL4 = (sum(pi.*(D4(D4>0).^4)./(128 * 1.002e-3)))^(-1)*length(D4(D4>0))./10^6; %MPa.s/m4
            RL = mean([RL1 RL2 RL3 RL4]);
            RW = RC - RL; %MPa.s/m4
            Rs = [RC RL RW];
            
        end
        
        function D = getD(obj,level)
            switch level
                case 1
                    yNode = 1;
                case 2
                    yNode = ceil(obj.Size(1)/3);
                case 3
                    yNode = floor(obj.Size(1)*2/3);
                case 4
                    D = sqrt(obj.Aend./pi).*2;
            end
            if level < 4
                x=yNode:obj.Size(1):(obj.Size(2)-1)*obj.Size(1)+yNode;
                D=zeros(1,obj.Size(2));
                for i=1:length(D)
                    temp=obj.getConduits(x(i));
                    if ~isempty(temp)
                        D(i) = temp.Diameter;
                    end
                end
            end
        end
        function [hCL,hCD,hLPP] = histograms(obj)
            %Returns conduit length, conduit diameter, and largest pit pore
            %diameter histograms
            CL = [obj.Conduits.Length].*10^3; %mm
            CD = [obj.Conduits.Diameter].*10^6; %micro meter
            LPP = [obj.ICConnections.Dpm].*10^9; %nm
            
            figure
            subplot(3,1,1)
            hCL=histogram(CL);
            %             hold on
            %             binMean=(hCL.BinEdges(2:end)+hCL.BinEdges(1:end-1))/2;
            %             f=fit(transpose(binMean),transpose(hCL.Values),'exp1');
            %             plot(f,binMean,transpose(hCL.Values))
            xlabel('Vessel Lengths (mm)')
            ylabel('Occurances')
            
            subplot(3,1,2)
            hCD=histogram(CD);
            xlabel('Vessel Diameters (\mum)')
            ylabel('Occurances')
            
            subplot(3,1,3)
            hLPP=histogram(LPP);
            xlabel('Largest Pit Pore Diameter (nm)')
            ylabel('Occurances')
        end
        
        function Embolized = get.Embolized(obj)
            if any(~[obj.Clusters.NonConducting])
                Embolized=false;
            else
                Embolized=true;
            end
        end
        
        function RepairAll(obj)
            for i=1:length(obj.Conduits)
                obj.Conduits(i).Repair;
                obj.Conduits(i).Integrate;
            end
        end
    end
    
end