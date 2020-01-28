classdef XylemNet < handle
    
    properties
        Conduits %Functional Conduits
        ConSorted
        ICConnections %Inter-Conduit Connections
        Kend %Vector containing last row of conductances
        %The following properties are useful for stitching
        Aend %End section vessel areas
        CBProw
        CBPcol
        pickedConx
        %graphs
        gCond
        gCav
        gBipNode
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
        
        function obj = XylemNet(row, column, depth,...
                Pe, NPe, Pc, Lce, Dc, Dc_cv, varargin)
            %Class constructor
            
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
            addRequired(reqInputs,'depth',@(x) validPosNum(x));
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
            parse(reqInputs, row, column, depth,...
                Pe, NPe, Pc, Lce, Dc, Dc_cv);
            
            obj.Size = [row column depth];
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
            obj.ICConnections = ...
                generateICCs(obj,Dp,Dm,k_ASP,lam_ASP,Fc,Fpf,Fap,Tm,Lp,...
                ASPcalcmethod);
            %Create graphs
            [obj.gCond, obj.gBipNode, obj.gCav] = createGraphs(obj);
            %updateConduits is a XylemNet method that deletes all
            %conduits that aren't part of any cluster and are therefore
            %not useful for the overall hydraulic pathway of the
            %segment
            obj.updateConduits();
            if isempty(obj.ICConnections)
                disp('No Network')
                return
            end
            
            %Sample conduit diameters from lognormal distribution
            Dstd=Dc_cv*Dc; %m
            Dm=log(Dc^2/sqrt(Dstd^2+Dc^2));
            Ds=sqrt(log(Dstd^2/(Dc^2)+1));
            Dcs=lognrnd(Dm,Ds,1,length(obj.Conduits));
            Dcs=sort(Dcs);
            tempLen=[obj.Conduits.Length];
            [~,ILen]=sort(tempLen);
            fn = DcMapping(length(obj.Conduits),...
                floor(length(obj.Conduits)/20));
            %Coordinate conduit diameter and length: n^th longest conduit
            %gets to be the n^th widest
            for k=1:length(obj.Conduits)
                obj.Conduits(k).updateDiameter(...
                    Dcs(fn(ILen==k)));
            end
            
            %Sort conduits according to their column positions for faster
            %searching. Here, the deleted conduits are removed.
            for i=1:length(obj.ConSorted)
                if isa(obj.ConSorted{i},'Conduit')
                    obj.ConSorted{i} =...
                        obj.ConSorted{i}(isvalid(obj.ConSorted{i}));
                end
            end
            %Figure out ID of connected conduits and statistics such as
            %number of connected conduits and number of inter-conduit
            %connections
            for i=1:length(obj.Conduits)
                obj.Conduits(i).addConConduits;
                obj.Conduits(i).updateConduit;
            end
            %Now area can be calculated and then conductance and similar
            %properties
            for i=1:length(obj.ICConnections)
                obj.ICConnections(i).updateMeanArea;
            end
            
            %Sample membrane stretching at air seeding (e) from a normal
            %distribution
            es = sort(abs(normrnd(e_mean,e_cv*e_mean,...
                1,length(obj.ICConnections))));
            %Distribute e according to ICC area and compute properties
            [~,I] = sort([obj.ICConnections.Am]);
            for i=1:length(obj.ICConnections)
                obj.ICConnections(i).computeKmASP(es(I==i));
            end
            %Figure out the minimum ASP with every connected conduit
            for i=1:length(obj.Conduits)
                obj.Conduits(i).addConConduitASP;
            end
            
            % update gCond Weights
            ceIdx = obj.gCond.Edges{:,"Type"}=="CE";
            iccIdx = obj.gCond.Edges{:,"Type"}=="ICC";
            
            ceObjs = obj.gCond.Edges{ceIdx,"CEObj"};
            iccObjs = obj.gCond.Edges{iccIdx,"ICCObj"};
            
            ceKs = [ceObjs.Kce];
            iccKs = [iccObjs.Km];
            
            obj.gCond.Edges{ceIdx,"Weight"} = ceKs';
            obj.gCond.Edges{iccIdx,"Weight"} = iccKs';
            
            % update gCav Weights
            for i = 1 : height(obj.gCav.Edges)
                endNodes = obj.gCav.Edges{i,"EndNodes"};
                cond1 = obj.gCav.Nodes{endNodes(1),"ConduitObj"};
                cond2 = obj.gCav.Nodes{endNodes(2),"ConduitObj"};
                obj.gCav.Edges{i, "Weight"} =...
                    cond1.ConConduitASP(cond1.ConConduits == cond2);
            end
            
            %Save organ border area and conductivity
            x = obj.Size(1):obj.Size(1):(obj.Size(2)*obj.Size(3))*obj.Size(1);
            obj.Kend=zeros(1,obj.Size(2)*obj.Size(3));
            obj.Aend=zeros(1,obj.Size(2)*obj.Size(3));
            for i=1:length(obj.Kend)
                temp = obj.getConduits(x(i));
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
            
            %Specify redundant conduits and corresponding nodes
            [~, condIdx] = redundancyCheck(obj.gCond, obj.gBipNode, obj.gCav);
            [obj.gCond, obj.gBipNode, obj.gCav] =...
                updateRedundancy(obj.gCond, obj.gBipNode, obj.gCav, condIdx);
        end
        
        %Network properties with embolization
        
        function updateConduits(obj)
            %This function removes non functional conduits, conduit
            %elements, and inter-conduit connections.
            condFunc = [obj.Conduits.Functional];
            iccFunc = [obj.ICConnections.Functional];
            %Delete non-functional Conduit elements
            rmCeIdx = find(~condFunc);
            for i = 1 : length(rmCeIdx)
               delete(obj.Conduits(rmCeIdx(i)).CEs); 
            end
            delete(obj.Conduits(~condFunc))
            delete(obj.ICConnections(~iccFunc))
            
            obj.Conduits = obj.Conduits(condFunc);
            obj.ICConnections = obj.ICConnections(iccFunc);
            
            for i = 1 : numel(obj.ConSorted)
               colDep = obj.ConSorted{i};
               if ~isempty(colDep)
                obj.ConSorted{i} = colDep(isvalid(colDep));
               end
            end
            
            %Delete non-functional ICCs
            for i = 1 : length(obj.Conduits)
               nonFunc = ~isvalid(obj.Conduits(i).ICConnections);
               obj.Conduits(i).ICConnections(nonFunc)= [];
            end
            
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
        
        function pg = visualize(obj, type)
            if ~exist('type', 'var') || type == "cond"
                pg = vizCond(obj.gCond);
            elseif type == "bip"
                pg = plot(obj.gBipNode);
            elseif type == "cav"       
                pg = vizCav(obj.gCav);
            end
        end
        
        function cons = getConduits(obj,Nodes)
            %Returns Conduits corresponding to the input Nodes
            % Only uses linear index
            cons=Conduit.empty();
            if size(Nodes,2) == 1 && size(Nodes,3) == 1
                 [I, J, K]= ind2sub(obj.Size,Nodes);
                 Nodes = [I; J; K]';
            end
            
            for j=1:size(Nodes,1)
                colDep = obj.ConSorted{1,Nodes(j,2),Nodes(j,3)};
                for i=1:length(colDep)
                        NN=colDep(i).Nodes(:,1);
                    if Nodes(j,1)>=NN(1) && Nodes(j,1)<=NN(end)
                        cons(j) = colDep(i);
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
                x = yNode : obj.Size(1) :...
                    (obj.Size(2)*obj.Size(3)-1)*obj.Size(1)+yNode;
                
                D=zeros(1,obj.Size(2)*obj.Size(3));
                
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
        
        function RepairAll(obj)
            for i=1:length(obj.Conduits)
                obj.Conduits(i).Repair;
                obj.Conduits(i).Integrate;
            end
        end
    end
    
end