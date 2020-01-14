%--------------------------------------------------------------------------
% FLWTFEM MATLAB SOLVER for FINITE ELEMENT ANALYSIS OF MULTILAYER PLATES
% BASED ON THE FULL LAYERWISE THEORY OF REDDY 
% Developed by: Miroslav Marjanovic, PhD Civil Eng.
%               Emilija Damnjanovic, MSc Civil Eng.
%               Belgrade, 2019.
%--------------------------------------------------------------------------
classdef OrthotropicLamina

    properties
        LaminaName = [];
        Material   = [];
        FiberAngle = 0;
        Thickness  = 0;
        
        % The Matrix T establishes the relationship between plate (XYZ) and
        % material (123) coordinate systems:
        % StressXYZ = T * Stress123;
        MatrixT    = [];
        
        % The Matrix R establishes the inverse relationship, between material
        % (123) and plate (XYZ) coordinate systems:
        % Stress123 = R * StressXYZ;
        MatrixR    = [];
        
        % The Matrix CBar establishes the relationship between stress and 
        % strain components in plate coordinates (XYZ):
        % StressXYZ = CBar * EpsilonXYZ
        MatrixCBar = zeros(6);
    end
        
    methods
        %Class Constructor
        function obj = OrthotropicLamina(name,mat,angle,thickness)
            obj.LaminaName = name;
            obj.Material   = mat;
            obj.FiberAngle = angle;
            obj.Thickness  = thickness;
            
            n = cosd(angle);
            m = sind(angle);
            
            obj.MatrixT = [ n^2  m^2  0   0   0   -2*m*n  ;
                            m^2  n^2  0   0   0    2*m*n  ;
                             0    0   1   0   0      0    ;
                             0    0   0   n   m      0    ;
                             0    0   0  -m   n      0    ;
                            m*n -m*n  0   0   0   n^2-m^2];
                         
            
            obj.MatrixR = [ n^2  m^2  0   0   0    2*m*n  ;
                            m^2  n^2  0   0   0   -2*m*n  ;
                             0    0   1   0   0      0    ;
                             0    0   0   n  -m      0    ;
                             0    0   0   m   n      0    ;
                           -m*n  m*n  0   0   0   n^2-m^2];              

            obj.MatrixCBar = obj.MatrixT * mat.MatrixCLocal * obj.MatrixT';
        end
    end
end