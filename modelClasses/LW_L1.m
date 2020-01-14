%--------------------------------------------------------------------------
% FLWTFEM MATLAB SOLVER for FINITE ELEMENT ANALYSIS OF MULTILAYER PLATES
% BASED ON THE FULL LAYERWISE THEORY OF REDDY 
% Developed by: Miroslav Marjanovic, PhD Civil Eng.
%               Emilija Damnjanovic, MSc Civil Eng.
%               Belgrade, 2019.
%--------------------------------------------------------------------------
classdef LW_L1
    
    properties
        ElementID;
        ElementType;
        Integration;
        Section;
        ElementNodes;
        DistributedLoad;
        NodalCoordinates;
        
        GaussPointsK;        GaussCoefsK;
        GaussPointsM;        GaussCoefsM;
        GaussPointsStress;   GaussCoefsStress;
    end
    
    methods
        
        % Class Constructor
        function obj = LW_L1(id,type,integration,section,nodes,load)
            if nargin~=0
                obj.ElementID       = id;
                obj.ElementType     = type;
                obj.Integration     = integration;
                obj.Section         = section;
                obj.ElementNodes    = nodes;
                obj.DistributedLoad = load;
                
                for i = 1:length(nodes)
                    obj.NodalCoordinates(i,1) = nodes(i).XCoord;
                    obj.NodalCoordinates(i,2) = nodes(i).YCoord;
                    obj.NodalCoordinates(i,3) = nodes(i).ZCoord;
                end
            end
        end
        
        % Derive Gauss Points and Coefficients
        function obj = deriveGauss(obj)
            
            % Triangular Elements = No Matter FULL or RED Integration
            % T3 or T6 = 3 GaussPoints for K, M and Stresses
            if strcmp(obj.ElementType,'T3') == 1 || strcmp(obj.ElementType,'T6') == 1
                [points3,coeffs3] = GaussCoordinates('Triangle', 3);
                
                obj.GaussPointsK = points3;
                obj.GaussCoefsK  = coeffs3;
                
                obj.GaussPointsM = points3;
                obj.GaussCoefsM  = coeffs3;
                
                obj.GaussPointsStress = points3;
                obj.GaussCoefsStress  = coeffs3;
            end

            
            % Q4  Element
            % FULL = 2x2 GaussPoints for K, M and also for Stresses
            % RED  = 2x2 GaussPoints for M and Stresses and 1 GaussPoint for K
            if strcmp(obj.ElementType,'Q4') == 1
                [points4,coeffs4] = GaussCoordinates('Quadrilateral', 4);
                [points1,coeffs1] = GaussCoordinates('Quadrilateral', 1);
                
                obj.GaussPointsK = points4;
                obj.GaussCoefsK  = coeffs4;
                
                obj.GaussPointsM = points4;
                obj.GaussCoefsM  = coeffs4;
                
                obj.GaussPointsStress = points4;
                obj.GaussCoefsStress  = coeffs4;
                    
                if strcmp(obj.Integration, 'RED') == 1
                    obj.GaussPointsK = points1;
                    obj.GaussCoefsK  = coeffs1;
                end 
            end
            
            % Q8  Element
            % FULL = 3x3 GaussPoints for K, M and 2x2 GaussPoints for Stresses
            % RED  = 3x3 GaussPoints for M and 2x2 GaussPoints for K and Stresses
            if strcmp(obj.ElementType,'Q8') == 1
                [points_9,coeffs_9] = GaussCoordinates('Quadrilateral', 9);
                [points4,coeffs4] = GaussCoordinates('Quadrilateral', 4);

                obj.GaussPointsK = points_9;
                obj.GaussCoefsK  = coeffs_9;
                    
                obj.GaussPointsM = points_9;
                obj.GaussCoefsM  = coeffs_9;
                    
                obj.GaussPointsStress = points4;
                obj.GaussCoefsStress  = coeffs4;
                    
                if strcmp(obj.Integration, 'RED') == 1
                    obj.GaussPointsK = points4;
                    obj.GaussCoefsK  = coeffs4;
                end
            end
            
        end
    end
end