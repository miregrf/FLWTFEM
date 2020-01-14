%--------------------------------------------------------------------------
% FLWTFEM MATLAB SOLVER for FINITE ELEMENT ANALYSIS OF MULTILAYER PLATES
% BASED ON THE FULL LAYERWISE THEORY OF REDDY 
% Developed by: Miroslav Marjanovic, PhD Civil Eng.
%               Emilija Damnjanovic, MSc Civil Eng.
%               Belgrade, 2019.
%--------------------------------------------------------------------------
classdef PostNode
    
    properties
        NodeID;
        XCoord;
        YCoord;
        ZCoord;
        
        DispX;
        DispY;
        DispZ;
    end
    
    
    methods
        % Class Constructor
        function obj = PostNode(id,x,y,z,varargin)
            if nargin ~=0
                obj.NodeID = id;
                obj.XCoord = x;
                obj.YCoord = y;
                obj.ZCoord = z;
            end
        end

    end
end