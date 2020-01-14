%--------------------------------------------------------------------------
% FLWTFEM MATLAB SOLVER for FINITE ELEMENT ANALYSIS OF MULTILAYER PLATES
% BASED ON THE FULL LAYERWISE THEORY OF REDDY 
% Developed by: Miroslav Marjanovic, PhD Civil Eng.
%               Emilija Damnjanovic, MSc Civil Eng.
%               Belgrade, 2019.
%--------------------------------------------------------------------------
classdef DOF
    properties
        ID;
        Free = 1;
    end
    
    methods
        % Class Constructor
        function obj = DOF(id)
            if nargin~=0
                obj.ID = id;
            end
        end
    end
end