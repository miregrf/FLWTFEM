%--------------------------------------------------------------------------
% FLWTFEM MATLAB SOLVER for FINITE ELEMENT ANALYSIS OF MULTILAYER PLATES
% BASED ON THE FULL LAYERWISE THEORY OF REDDY 
% Developed by: Miroslav Marjanovic, PhD Civil Eng.
%               Emilija Damnjanovic, MSc Civil Eng.
%               Belgrade, 2019.
%--------------------------------------------------------------------------
classdef Node
    
    properties
        NodeID;
        XCoord;
        YCoord;
        ZCoord;

        % Assigned Nodal Forces
        ForceX_top = 0;
        ForceY_top = 0;
        ForceZ_top = 0;
        ForceX_bottom = 0;
        ForceY_bottom = 0;
        ForceZ_bottom = 0;
        
        % List of Nodal DOFs
        NodalDOFs = [];
    end
    
    methods
        % Class Constructor
        function obj = Node(id,coordinates,varargin)
            if nargin ~=0
                obj.NodeID = id;
                obj.XCoord = coordinates(obj.NodeID,1);
                obj.YCoord = coordinates(obj.NodeID,2);
                obj.ZCoord = coordinates(obj.NodeID,3);
            end
        end
        
        % Assignment of Nodal Forces in Nodes
        function obj = assignNodalForces(obj, force)
            obj.ForceX_top = force(1);
            obj.ForceY_top = force(2);
            obj.ForceZ_top = force(3);
            obj.ForceX_bottom = force(4);
            obj.ForceY_bottom = force(5);
            obj.ForceZ_bottom = force(6);
        end
        
        % Assignment of DOFs in Nodes
        function obj = assignDOFs(obj, NDOF, DOFS)
            obj.NodalDOFs = DOFS ( (obj.NodeID-1)*NDOF+1 : obj.NodeID*NDOF );
        end
    end
end