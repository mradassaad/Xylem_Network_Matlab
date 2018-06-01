classdef ICC < handle
    %EndWall supports different modeling schemes:
    %1 - Takes as inputs Dp,Dp_cv,Fc,Fpf
    properties
        Dpm %maximum pit pore diameter meters (m)
        Nodes
        OINodes
        ConConduits
        
        Npore %number of pores in connection
        Npit %number of pits in connection
        Am %inter-vessel area occupied by pits
        Dp %mean pore diameter (m)
        Dp_cv %coefficient of variation of pore diameters
        De %effective hydraulic pore diameters (m)
        Dm %average diameter of pit membranes
        Km %Hydraulic conductance %m3.MPa-1.s-1
        Tm %pit membrane thicness (m)
        A %pit air-Seeding scale parameter %MPa
        B %Shape parameter
        Pb %Pressure at which membrane touches chamber walls
        ASP %endwall air-seeding pressure %MPa
        Fc %fraction of vessel surface area in contact with other vessels
        Fpf %fraction of contact area occupied by pits (pit-field)
        Fap %fraction of aperture area to membrane area
        Lp %Pit chamber depth (m)
        emax %Maximum membrane strain at air seeding

        Embolized %Specifies if embolized
        NonConducting %Specifies if isolated
        ASPcalcmethod %'Pore' or 'Stretching' 
    end
    
    properties (Constant)
        mu = 1.002e-3; %Pa.sd
        gamma = 0.072; %Pa
        E = 400; % Young's modulus of elasticity (MPa)
        v = 0.3; %Poisson ratio
    end
    
    methods
        
        function obj = ICC(Nodes,ConConduits,dimensions,...
                ASPcalcmethod,varargin)
            if nargin<=0
                obj.Dp=0;
                obj.Dp_cv = 0;
                obj.Dpm=0;
                obj.De=0;
                obj.Nodes=0;
                obj.ConConduits = Conduit.empty;
                obj.Embolized=false;
                obj.NonConducting=false;
                obj.Fc = 0;
                obj.Fpf = 0;
                obj.Npore = 0;
                obj.Npit = 0;
                obj.A = 0;
                obj.B = 0;
                obj.Dm = 0;
                obj.Tm = 0;
            else
                obj.Nodes=Nodes;
                obj.ConConduits = ConConduits; %The connected conduits
                obj.OINodes=obj.OneIndexNodes(dimensions);
                obj.ASPcalcmethod = ASPcalcmethod;
                obj.Embolized = false;
                obj.NonConducting = false;

                obj.Dp = varargin{1};
                obj.Dp_cv = 0;
                obj.Fc = varargin{5};
                obj.Fpf = varargin{6};
                obj.Fap = varargin{7};
                obj.A = varargin{3};
                obj.B = varargin{4};
                obj.Dm = mean(varargin{2});
                obj.Tm = mean(varargin{8});
                obj.Lp = mean(varargin{9});
                obj.Dpm = 0;
                obj.De = obj.Dp;
                obj.Am = 0;
                obj.Km = 0;
                obj.Npore = 0;
                obj.Npit = 0;
            end
            
        end
        function updateMeanArea(obj)
            obj.Am = (1/2)*...
                (obj.ConConduits(1).Ac/obj.ConConduits(1).ICCCount...
                + obj.ConConduits(2).Ac/obj.ConConduits(2).ICCCount)...
                *obj.Fc*obj.Fpf;
                %total area occupied by a pore with microfibril strand
                tf = 30e-9; %m - microfibril strand thickness (sperry hacke)
                Apore = pi*(obj.Dp+tf)^2/4; %area of single pore (m2)
                obj.Npore = floor(obj.Am/(Apore));

                Apit = obj.Dm^2; %area of single pit (m2)
                obj.Npit = floor(obj.Am/Apit);
                
        end
        
        function computeKmASP(obj,emax)
            obj.emax = emax;
            %Assign the following outous depending on which method is
            %used for ASP calculation
            [ASPout,Pbout] = deflectionBPP(obj.Dm,obj.Tm,obj.v,obj.E,...
                obj.emax,obj.Lp,obj.Fap);
            obj.Pb = Pbout;
            if isequal(obj.ASPcalcmethod,'Pore')
                obj.ASP = wblrnd(obj.A/obj.Npit^(1/obj.B),obj.B);
            elseif isequal(obj.ASPcalcmethod,'Thickness')
                obj.ASP = ASPout;
            end
            obj.Dpm = 4*obj.gamma/obj.ASP * 1e-6; %m
            
            Kmtemp = (obj.De^3/(24*obj.mu))*...
                (1+16*obj.Tm/(3*pi*obj.De))^(-1)*obj.Npore; %m3.Pa-1.s-1
            obj.Km = Kmtemp*1e6; %m3.MPa-1.s-1
        end
        
        function OINodes = OneIndexNodes(obj,dimensions)
            %Convert conduit and pit matrices to equivalent positions in the new
            %numbering system
%             OINodes=dimensions(1)*(obj.Nodes(:,2)-1)+obj.Nodes(:,1);
            OINodes=sub2ind(dimensions,obj.Nodes(:,1),obj.Nodes(:,2));
        end
        
        function OINodesShift = ShiftedOneIndexNodes(obj,dimensions,clusterColumns)
            %Convert conduit and pit matrices to equivalent positions in the new
            %numbering system inside clusters
            OINodesShift=dimensions(1)*((obj.Nodes(:,2)-clusterColumns(1)+1)-1)+obj.Nodes(:,1);
        end
 
        function Embolize(obj)
            obj.Embolized=true;
            obj.NonConducting=true;
        end
        
        function Isolate(obj)
            obj.NonConducting=true;
        end
        
        function Repair(obj)
            if all(~[obj.ConConduits.Embolized])
                obj.Embolized=false;
            end
        end
        
        function Integrate(obj)
            if all(~[obj.ConConduits.NonConducting])
                obj.NonConducting=false;
            end
        end
        
    end
end