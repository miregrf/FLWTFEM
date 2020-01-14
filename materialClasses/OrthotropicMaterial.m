%--------------------------------------------------------------------------
% FLWTFEM MATLAB SOLVER for FINITE ELEMENT ANALYSIS OF MULTILAYER PLATES
% BASED ON THE FULL LAYERWISE THEORY OF REDDY 
% Developed by: Miroslav Marjanovic, PhD Civil Eng.
%               Emilija Damnjanovic, MSc Civil Eng.
%               Belgrade, 2019.
%--------------------------------------------------------------------------
classdef OrthotropicMaterial
    
    properties
        MaterialName = [];
        
        E1 = 0; 
        E2 = 0; 
        E3 = 0;
        nu12 = 0; 
        nu13 = 0; 
        nu23 = 0;
        G12 = 0; 
        G13 = 0; 
        G23 = 0;
        Rho = 0;
        
        nu21 = 0; 
        nu31 = 0; 
        nu32 = 0;
        
        % MatrixCLocal establishes the relationship between stress and 
        % strain components in material coordinates (123):
        % Stress123 = MatrixCLocal * Epsilon123
        MatrixCLocal;
    end
    
    methods
        % Class Constructor
        function obj = OrthotropicMaterial(name, E1,E2,E3,nu12,nu13,nu23,G12,G13,G23, Rho)
            obj.MaterialName = name;
            obj.E1 = E1;
            obj.E2 = E2;
            obj.E3 = E3;
            obj.nu12 = nu12;
            obj.nu13 = nu13;
            obj.nu23 = nu23;
            obj.G12 = G12;
            obj.G13 = G13;
            obj.G23 = G23;
            obj.Rho = Rho;
            
            obj.nu21 = E2*nu12/E1;
            obj.nu31 = E3*nu13/E1;
            obj.nu32 = E3*nu23/E2;
            
            delta = (1 - nu12*obj.nu21 - nu23*obj.nu32 - nu13*obj.nu31 - 2*obj.nu21*obj.nu32*nu13) / E1 / E2 / E3;
            C11   = (1-nu23*obj.nu32) / E2/E3/delta;
            C12   = (obj.nu21 + obj.nu31*nu23) / E2/E3/delta;
            C13   = (obj.nu31 + obj.nu21*obj.nu32) / E2/E3/delta;
            C22   = (1 - nu13*obj.nu31) / E1/E3/delta;
            C23   = (obj.nu32 + nu12*obj.nu31) / E1/E3/delta;
            C33   = (1 - nu12*obj.nu21) / E1/E2/delta;
            
            C44 = G23;
            C55 = G13;
            C66 = G12;
            
            obj.MatrixCLocal = [C11 C12 C13   0   0   0;
                                C12 C22 C23   0   0   0;
                                C13 C23 C33   0   0   0;
                                 0   0   0   C44  0   0;
                                 0   0   0   0   C55  0;
                                 0   0   0   0    0  C66];
        end
    end
end      