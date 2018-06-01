classdef CondElem < handle
    properties
        Diameter
        Length
        Nodes %Nodes connected where [a1 b1;a2 b2]
        Embolised
        Kce %Conductivity of conduit element in m^3/MPa.s
    end

    properties (Constant)
        mu = 1.002e-3; %Pa.s
    end
    
    methods
        function obj = CondElem(len,nodes)
            if nargin<=0
                obj.Diameter=0;
                obj.Length=0;
                obj.Embolised=false;
                obj.Nodes=zeros(2,2);
                obj.Kce=0;
            else
                %make sure diameter and length are positive numbers
                if len<=0
                    error('Diameter and Length of conduit should be positive numbers')
                end
                
                obj.Diameter=0; %m
                obj.Length=len; %m
                %Make sure nodes are different
                if nodes(1,1)==nodes(2,1) && nodes(1,2)==nodes(2,2)
                    error('Nodes should be different')
                end
                %Make sure that inputted nodes are consecutive nodes either
                %diagonal or axial
                if ~(abs(nodes(1,1)-nodes(2,1))<=1 && abs(nodes(1,2)-nodes(2,2))<=1)
                    error('Nodes in conduit element must be consecutive')
                end
                %Make sure conduit cannot be horizontal
                if nodes(1,1)==nodes(2,1) && abs(nodes(1,2)-nodes(2,2))==1
                    error('Conduit element must not be horizontal')
                end
                obj.Nodes=nodes;
                obj.Embolised=false;
                obj.Kce = 0;
            end
        end 
        function updateDiameter(obj,dia)
            obj.Diameter=dia; %m
            res= 128*obj.mu*obj.Length/(pi*obj.Diameter^4);
            obj.Kce = 1/res; %m^3/Pa.s
            obj.Kce = obj.Kce*10^6; %m^3/MPa.s
        end
    end
end