%--------------------------------------------------------------------------
% FLWTFEM MATLAB SOLVER for FINITE ELEMENT ANALYSIS OF MULTILAYER PLATES
% BASED ON THE FULL LAYERWISE THEORY OF REDDY 
% Developed by: Miroslav Marjanovic, PhD Civil Eng.
%               Emilija Damnjanovic, MSc Civil Eng.
%               Belgrade, 2019.
%--------------------------------------------------------------------------
classdef PostElement
    
    properties
        ElementID;
        ElementType;
        ElementNodes;
        
        StressesXYZ;
        Stresses123;
    end
    
    methods
        % Class Constructor
        function obj = PostElement(id, postnodes)
            if nargin~=0
                obj.ElementID = id;
                obj.ElementNodes = postnodes;
                
                if length(postnodes) == 6
                    obj.ElementType = 'Prism 6';
                elseif length(postnodes) == 8
                    obj.ElementType = 'Hexahedra 8';
                end
            end
        end
    end
end