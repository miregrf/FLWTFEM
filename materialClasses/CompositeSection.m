%--------------------------------------------------------------------------
% FLWTFEM MATLAB SOLVER for FINITE ELEMENT ANALYSIS OF MULTILAYER PLATES
% BASED ON THE FULL LAYERWISE THEORY OF REDDY 
% Developed by: Miroslav Marjanovic, PhD Civil Eng.
%               Emilija Damnjanovic, MSc Civil Eng.
%               Belgrade, 2019.
%--------------------------------------------------------------------------
classdef CompositeSection
    
    properties
        SectionName;
        Layers;
        MidLayer;
        MatrixL;

        Inertia_IJ;
        
        MatrixA1;
        MatrixA2;
    end
    
    methods
        % Class Constructor
        function obj = CompositeSection(name, laminas)
            obj.SectionName = name;
            obj.Layers = laminas;
            
            zz = zeros(length(obj.Layers) + 1, 1);
            hh = 0;
            for i = 2:length(zz)
                hh = hh + obj.Layers(i-1).Thickness;
                zz(i) = hh;
            end
            
            hh = zz(1);
            for i = 1:length(obj.Layers)
                hh = hh + obj.Layers(i).Thickness;
                if hh > -1e-9 && hh < 1e-9
                    obj.MidLayer = i+1;
                end
            end
            
            %MatrixL stores geometric properties of the Layers within the
            %considered section, for the calculation of interlaminar
            %stresses
            N = size(obj.Layers,2); %Number of Material Layers
            L = zeros(3*N, 3*N);
            
            L(1, 3)       = 1;
            L(2, 3*N-2)   = (obj.Layers(N).Thickness)^2;
            L(2, 3*N-1)   = obj.Layers(N).Thickness;
            L(2, 3*N)     = 1;
            
            for I = 3:N+1
                J = I-2;
                L(I,3*J-2)    = (obj.Layers(J).Thickness)^2;
                L(I,3*J-2+1)  = obj.Layers(J).Thickness;
                L(I,3*J-2+2)  = 1;
                L(I,3*J-2+5)  = -1; 
            end
            
            for I = N+2:2*N+1
                J = I-N-1;
                L(I,3*J-2)    = (obj.Layers(J).Thickness)^3/3;
                L(I,3*J-2+1)  = (obj.Layers(J).Thickness)^2/2;
                L(I,3*J-2+2)  = obj.Layers(J).Thickness;
            end
            
            for I = 2*N+2:3*N
                J = I-2*N-1;
                L(I,3*J-2)    = 2*obj.Layers(J).Thickness;
                L(I,3*J-2+1)  = 1;
                L(I,3*J-2+4)  = -1;            
                
            end              
            obj.MatrixL = L;
        end
        
        
        % Calculation of Characteristic Coefficient Matrices for 3D Stress Case
        function obj = calcCoeffMatrices(obj)
            N = length(obj.Layers); %Number of Material Layers
            Apq       = zeros(6,6, N+1, N+1);
            ApqTilde  = zeros(6,6, N+1, N+1);
            ApqBar    = zeros(6,6, N+1, N+1);
            ApqBarBar = zeros(6,6, N+1, N+1);
            
            Apq(:,:, 1,   1)   = obj.Layers(1).MatrixCBar * obj.Layers(1).Thickness / 3;
            Apq(:,:, 1,   2)   = obj.Layers(1).MatrixCBar * obj.Layers(1).Thickness / 6;
            Apq(:,:, N+1, N)   = obj.Layers(N).MatrixCBar * obj.Layers(N).Thickness / 6;
            Apq(:,:, N+1, N+1) = obj.Layers(N).MatrixCBar * obj.Layers(N).Thickness / 3;
            
            ApqTilde(:,:, 1,   1)   = -obj.Layers(1).MatrixCBar / 2;
            ApqTilde(:,:, 1,   2)   =  obj.Layers(1).MatrixCBar / 2;
            ApqTilde(:,:, N+1, N)   = -obj.Layers(N).MatrixCBar / 2;
            ApqTilde(:,:, N+1, N+1) =  obj.Layers(N).MatrixCBar / 2;
            
            ApqBar(:,:, 1,   1)   =  obj.Layers(1).MatrixCBar / obj.Layers(1).Thickness;
            ApqBar(:,:, 1,   2)   = -obj.Layers(1).MatrixCBar / obj.Layers(1).Thickness;
            ApqBar(:,:, N+1, N)   = -obj.Layers(N).MatrixCBar / obj.Layers(N).Thickness;
            ApqBar(:,:, N+1, N+1) =  obj.Layers(N).MatrixCBar / obj.Layers(N).Thickness;
            
            ApqBarBar(:,:, 1,   1)   = -obj.Layers(1).MatrixCBar / 2;
            ApqBarBar(:,:, 1,   2)   = -obj.Layers(1).MatrixCBar / 2;
            ApqBarBar(:,:, N+1, N  ) = obj.Layers(N).MatrixCBar / 2;
            ApqBarBar(:,:, N+1, N+1) = obj.Layers(N).MatrixCBar / 2;

            for I = 2:N
                Apq(:,:, I, I-1) = obj.Layers(I-1).MatrixCBar * obj.Layers(I-1).Thickness / 6;
                Apq(:,:, I, I)   = obj.Layers(I-1).MatrixCBar * obj.Layers(I-1).Thickness / 3 + obj.Layers(I).MatrixCBar * obj.Layers(I).Thickness / 3;
                Apq(:,:, I, I+1) = obj.Layers(I).MatrixCBar * obj.Layers(I).Thickness / 6;
                
                ApqTilde(:,:, I, I-1) = -obj.Layers(I-1).MatrixCBar / 2;
                ApqTilde(:,:, I, I)   =  obj.Layers(I-1).MatrixCBar / 2 - obj.Layers(I).MatrixCBar / 2;
                ApqTilde(:,:, I, I+1) =  obj.Layers(I).MatrixCBar / 2;
                
                ApqBar(:,:, I, I-1) = -obj.Layers(I-1).MatrixCBar / obj.Layers(I-1).Thickness;
                ApqBar(:,:, I, I)   =  obj.Layers(I-1).MatrixCBar/obj.Layers(I-1).Thickness   +   obj.Layers(I).MatrixCBar/obj.Layers(I).Thickness;
                ApqBar(:,:, I, I+1) = -obj.Layers(I).MatrixCBar / obj.Layers(I).Thickness;
                
                ApqBarBar(:,:, I, I-1) =  obj.Layers(I-1).MatrixCBar / 2;
                ApqBarBar(:,:, I, I)   =  obj.Layers(I-1).MatrixCBar / 2 - obj.Layers(I).MatrixCBar / 2;
                ApqBarBar(:,:, I, I+1) = -obj.Layers(I).MatrixCBar / 2;
            end
            
            A1 = zeros(4,4,N+1,N+1);
            A2 = zeros(4,4,N+1,N+1);
            for I = 1:N+1
                for J = 1:N+1
                    A1_11 = Apq([1 2 6], [1 2 6], I, J);
                    A1_12 = ApqTilde([1 2 6], 3, I, J);
                    A1_21 = ApqBarBar(3, [1 2 6], I, J);
                    A1_22 = ApqBar(3, 3, I, J);
                    A1(:,:,I,J) = [A1_11 A1_12; A1_21 A1_22];
                    
                    A2_11 = ApqBar([5 4], [5 4], I, J);
                    A2_12 = ApqBarBar([5 4], [5 4], I, J);
                    A2_21 = ApqTilde([5 4], [5 4], I, J);
                    A2_22 = Apq([5 4], [5 4], I, J);
                    A2(:,:,I,J) = [A2_11 A2_12; A2_21 A2_22];
                end
            end
            
            obj.MatrixA1 = A1;
            obj.MatrixA2 = A2;
        end

        
        % Calculation of Characteristic Inertia Matrices for 3D Stress Case
        function obj = calcInertiaMatrices(obj)
            N = length(obj.Layers); %Number of Material Layers
            InertiaMatrix = zeros(N+1, N+1);
            
            InertiaMatrix(1,1)       = obj.Layers(1).Material.Rho/3 * obj.Layers(1).Thickness;
            InertiaMatrix(1,2)       = obj.Layers(1).Material.Rho/6 * obj.Layers(1).Thickness;
            InertiaMatrix(N+1,N)     = obj.Layers(N).Material.Rho/6 * obj.Layers(N).Thickness;
            InertiaMatrix(N+1,N+1)   = obj.Layers(N).Material.Rho/3 * obj.Layers(N).Thickness;
            for I = 2:N
                InertiaMatrix(I,I-1) = obj.Layers(I-1).Material.Rho/6 * obj.Layers(I-1).Thickness;
                InertiaMatrix(I,I)   = (obj.Layers(I-1).Material.Rho  * obj.Layers(I-1).Thickness + obj.Layers(I).Material.Rho * obj.Layers(I).Thickness) / 3;
                InertiaMatrix(I,I+1) = obj.Layers(I).Material.Rho/6   * obj.Layers(I).Thickness;
            end
            obj.Inertia_IJ = InertiaMatrix;
        end
        
    end
end              