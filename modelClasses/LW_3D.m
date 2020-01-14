%--------------------------------------------------------------------------
% FLWTFEM MATLAB SOLVER for FINITE ELEMENT ANALYSIS OF MULTILAYER PLATES
% BASED ON THE FULL LAYERWISE THEORY OF REDDY 
% Developed by: Miroslav Marjanovic, PhD Civil Eng.
%               Emilija Damnjanovic, MSc Civil Eng.
%               Belgrade, 2019.
%--------------------------------------------------------------------------
classdef LW_3D < LW_L1
    
    properties
        MatrixKL;
        MatrixM;
        VectorQ;
        
        IndexDOF;
        Displacements;
        
        StressesXYZ;
        Stresses123
    end
    
    methods
        
        % Class Constructor
        function obj = LW_3D(id,type,integration,section,nodes,load)
            obj = obj@LW_L1(id,type,integration,section,nodes,load);
            if nargin~=0
            end
        end
        
        
        % Calculate Element Stiffness Matrix
        function obj = calcMatrixKL(obj)
            
            nnodes  = length(obj.ElementNodes);
            nodedoflayer = 3*nnodes;
            nlayers = length(obj.Section.Layers);
            nq      = length(obj.GaussCoefsK);
            K       = zeros ( 3 * nnodes * (nlayers+1) );
            
            for gp = 1:nq
                xi = obj.GaussPointsK(gp,1);
                eta = obj.GaussPointsK(gp,2);
                
                [shape,NatDev,~] = Shape(obj.ElementType, xi, eta);
                
                J11 = NatDev(:,1)' * obj.NodalCoordinates(:,1);
                J12 = NatDev(:,1)' * obj.NodalCoordinates(:,2);
                J21 = NatDev(:,2)' * obj.NodalCoordinates(:,1);
                J22 = NatDev(:,2)' * obj.NodalCoordinates(:,2);
                detJ = J11*J22 - J12*J21;
                invJ = 1/detJ * [J22 -J21; -J12 J11];
                XYDev = NatDev * invJ;
                
                matrixB(1:3, :)    = kron(XYDev(:,1)', [1 0 0; 0 0 0; 0 1 0]) + kron(XYDev(:,2)', [0 0 0; 0 1 0;1 0 0]);
                matrixB(4,   :)    = kron(shape, [0 0 1]);
                
                matrixBbar(1:2, :) = kron(shape,  [1 0 0; 0 1 0]);
                matrixBbar(3:4, :) = kron(XYDev', [0 0 1]);
                
                for I = 1:nlayers + 1
                    for J = 1:nlayers + 1
                        if abs(I-J) < 2
                            
                            K ( (I-1)*nodedoflayer+1:I*nodedoflayer , (J-1)*nodedoflayer+1:J*nodedoflayer )  =  K ( (I-1)*nodedoflayer+1:I*nodedoflayer , (J-1)*nodedoflayer+1:J*nodedoflayer ) + (matrixB' * obj.Section.MatrixA1(:,:,I,J) * matrixB   +   matrixBbar' * obj.Section.MatrixA2(:,:,I,J) * matrixBbar) * detJ * obj.GaussCoefsK(gp);
                        
                        end
                    end
                end
            end
            obj.MatrixKL = K;
        end
        
        
        % Calculate Element Mass Matrix
        function obj = calcMatrixM(obj)
            nnodes  = length(obj.ElementNodes);
            nodedoflayer = 3*nnodes;
            nlayers = length(obj.Section.Layers);
            nq      = length(obj.GaussCoefsK);
            M       = zeros ( 3 * nnodes * (nlayers+1) );
            
            for gp = 1:nq
                xi = obj.GaussPointsM(gp,1);
                eta = obj.GaussPointsM(gp,2);
                
                [shape,NatDev,~] = Shape(obj.ElementType, xi, eta);
                
                J11 = NatDev(:,1)' * obj.NodalCoordinates(:,1);
                J12 = NatDev(:,1)' * obj.NodalCoordinates(:,2);
                J21 = NatDev(:,2)' * obj.NodalCoordinates(:,1);
                J22 = NatDev(:,2)' * obj.NodalCoordinates(:,2);
                detJ = J11*J22 - J12*J21;
                
                matrixN  = kron(shape, eye(3));
                PRODUCT = matrixN' * matrixN * detJ;
                
                for I = 1:nlayers + 1
                    for J = 1:nlayers + 1
                        if abs(I-J) < 2
                            
                            M ( (I-1)*nodedoflayer+1:I*nodedoflayer , (J-1)*nodedoflayer+1:J*nodedoflayer )  =  M ( (I-1)*nodedoflayer+1:I*nodedoflayer , (J-1)*nodedoflayer+1:J*nodedoflayer )  +  obj.Section.Inertia_IJ(I,J) * PRODUCT * obj.GaussCoefsM(gp);
                        
                        end
                    end
                end  
                
            end
            obj.MatrixM = M;
        end

        
        % Calculate Element Force Vector
        function obj = calcVectorQ(obj)
            nnodes = length(obj.ElementNodes);
            nlayers = length(obj.Section.Layers);
            
            Q = zeros(nnodes * 3 * (nlayers + 1), 1);
            
            for gp = 1:length(obj.GaussCoefsM)
                xi = obj.GaussPointsM(gp,1);
                eta = obj.GaussPointsM(gp,2);
                
                [shape,NatDev,~] = Shape(obj.ElementType, xi, eta);
                
                J11 = NatDev(:,1)' * obj.NodalCoordinates(:,1);
                J12 = NatDev(:,1)' * obj.NodalCoordinates(:,2);
                J21 = NatDev(:,2)' * obj.NodalCoordinates(:,1);
                J22 = NatDev(:,2)' * obj.NodalCoordinates(:,2);
                detJ = J11*J22 - J12*J21;
                
                PRODUCT1 = ( kron(shape, eye(3)) )'  *  obj.DistributedLoad(4:6)'  *  detJ;
                PRODUCT2 = ( kron(shape, eye(3)) )'  *  obj.DistributedLoad(1:3)'  *  detJ;
                
                Q(1 : 3*nnodes , 1)               = Q(1:3*nnodes,1)               + PRODUCT1 * obj.GaussCoefsM(gp);
                Q(3*nnodes*nlayers + 1 : end , 1) = Q(3*nnodes*nlayers + 1:end,1) + PRODUCT2 * obj.GaussCoefsM(gp);
            end
            obj.VectorQ = Q;
        end
        
        
        % Derive Element Code Numbers
        function obj = deriveIndexDOF(obj)
            nnodes = length(obj.ElementNodes);
            nlayers = length(obj.Section.Layers);
            ndof = 3 * (nlayers+1);
            
            obj.IndexDOF = zeros( nnodes * ndof , 1 );
            for layer = 1:nlayers + 1
                for node = 1:nnodes
                    obj.IndexDOF((layer-1)*3*nnodes + 3*node - 2) = (obj.ElementNodes(node).NodeID-1)*ndof + 3*layer - 2;
                    obj.IndexDOF((layer-1)*3*nnodes + 3*node - 1) = (obj.ElementNodes(node).NodeID-1)*ndof + 3*layer - 1;
                    obj.IndexDOF((layer-1)*3*nnodes + 3*node    ) = (obj.ElementNodes(node).NodeID-1)*ndof + 3*layer;
                end
            end
        end
        
        
        % Calculate Element Stresses
        function obj = CalcStresses(obj, sublayers)
            ngaus   = size(obj.GaussPointsStress,1);
            nlayers = length(obj.Section.Layers);
            nnodes  = length(obj.ElementNodes);
            
            %Rows    = Nodes
            %Columns = X Y
            NodalCoords = obj.NodalCoordinates(:,1:2); 
            
            % PreliminaryStressMatrixGaussXYZ
            % Rows    = Six Componental Stresses SigmaXX SigmaYY SigmaZZ TauYZ TauXZ TauXY
            % Columns = Two Interfaces (Bottom and Top) of the Current Elem
            % Third   = Layer
            % Fourth  = Gauss Points
            PreliminaryStressMatrixGaussXYZ = zeros(6,2,nlayers,ngaus);          
            
            % derStressMatrixGaussXYZ
            % Rows    = Three Derivatives dTauXZ/dX dTauYZ/dY dSigmaZZ/dZ
            % Columns = Two Interfaces (Bottom and Top) of the Current Elem
            % Third   = Layer
            % Fourth  = Gauss Points
            derStressMatrixGaussXYZ = zeros(3,2,nlayers,ngaus);  
            
            % TauVectorGaussXYZ
            % Rows    = 3*n Values associated with Material Layers
            % Columns = TauXZ TauYZ SigmaZZ
            % Third   = Gauss Points
            TauVectorGaussXYZ = zeros(3*nlayers,3,ngaus);
            
            % Assignment of Traction Boundary Conditions on the Top and Bottom Surfaces
            TauVectorGaussXYZ(1,1,:) =  obj.DistributedLoad(4); %TauXZ - Bottom
            TauVectorGaussXYZ(2,1,:) =  obj.DistributedLoad(1); %TauXZ - Top
            TauVectorGaussXYZ(1,2,:) =  obj.DistributedLoad(5); %TauYZ - Bottom
            TauVectorGaussXYZ(2,2,:) =  obj.DistributedLoad(2); %TauYZ - Top
            TauVectorGaussXYZ(1,3,:) = -obj.DistributedLoad(6); %SigmaZZ - Bottom
            TauVectorGaussXYZ(2,3,:) =  obj.DistributedLoad(3); %SigmaZZ - Top
            
            % VectorCGaussXYZ
            % Rows    = 3*n Coefficients for Quadratic Interpolation of Stresses through a Single Material Layer in a Considered Gauss Point
            % Columns = TauXZ TauYZ SigmaZZ
            % Third   = Gauss Points
            VectorCGaussXYZ = zeros(3*nlayers, 3, ngaus);

            % StressMatrixGaussXYZ
            % Rows    = Six Componental Stresses SigmaXX SigmaYY SigmaZZ TauYZ TauXZ TauXY
            % Columns = Two Interfaces (Bottom and Top) of the Current Elem
            % Third   = Layer
            % Fourth  = Gauss Points
            StressMatrixGaussXYZ = zeros(6, 2, nlayers, ngaus);
            
            % StressesXYZ
            % Rows    = Six Componental Stresses SigmaXX SigmaYY SigmaZZ TauYZ TauXZ TauXY - Final Values
            % Columns = Interfaces (sublayers+1) of the Current Elem
            % Third   = Layer
            % Fourth  = Gauss Points
            obj.StressesXYZ = zeros(6, sublayers+1, nlayers, ngaus);
            obj.Stresses123 = zeros(6, sublayers+1, nlayers, ngaus);
            
                
            % Loop Through all Numerical Layers
            for layer = 1:nlayers
                
                % Extract MatrixCBar and Thickness of the Current Numerical Layer of the Finite Element
                CBar = obj.Section.Layers(layer).MatrixCBar;
                h    = obj.Section.Layers(layer).Thickness;
                
                %Extract Nodal Displacements at the Bottom and the Top of the Current Numerical Layer of the Finite Element
                Disp_Bottom = obj.Displacements( (layer-1)*3*nnodes+1 :   layer  *3*nnodes);
                Disp_Top    = obj.Displacements(  layer   *3*nnodes+1 : (layer+1)*3*nnodes);
                
                % Loop Gauss Points
                for gp = 1:ngaus
                    
                    xi = obj.GaussPointsStress(gp,1);
                    eta = obj.GaussPointsStress(gp,2);
                    
                    [shape,NatDev,SecDev] = Shape(obj.ElementType, xi, eta);
                    
                    J11 = NatDev(:,1)' * obj.NodalCoordinates(:,1);
                    J12 = NatDev(:,1)' * obj.NodalCoordinates(:,2);
                    J21 = NatDev(:,2)' * obj.NodalCoordinates(:,1);
                    J22 = NatDev(:,2)' * obj.NodalCoordinates(:,2);
                    detJ = J11*J22 - J12*J21;
                    invJ = 1/detJ * [J22 -J21; -J12 J11];
                    XYDev = NatDev * invJ;
                
                    J1 = [J11^2   J12^2   2*J11*J12;
                          J21^2   J22^2   2*J21*J22;
                          J11*J21 J12*J22 J21*J12 + J11*J22];
                    
                    J2 = SecDev' * NodalCoords;
                    
                    MatrixB1    = zeros(3, 3*nnodes);
                    MatrixB2    = kron(shape, [0 0 1]);
                    MatrixBbar1 = kron(shape, [1 0 0;0 1 0]);
                    MatrixBbar2 = zeros(2, 3*nnodes);
                    
                    MatrixH1    = zeros(3, 3*nnodes);
                    MatrixH2    = zeros(3, 3*nnodes);
                    MatrixG1    = zeros(1, 3*nnodes);
                    MatrixG2    = zeros(1, 3*nnodes);
                    MatrixM1    = zeros(2, 3*nnodes);
                    MatrixM2    = zeros(2, 3*nnodes);
                    MatrixN1    = zeros(1, 3*nnodes);
                    MatrixN2    = zeros(1, 3*nnodes);
                    
                    for i = 1:nnodes
                        MatrixB1   (:, 3*i-2:3*i) = [XYDev(i,1) 0 0; 0 XYDev(i,2) 0; XYDev(i,2) XYDev(i,1) 0];
                        MatrixBbar2(:, 3*i-2:3*i) = [0 0  XYDev(i,1); 0 0  XYDev(i,2)];
                        
                        dNd2 = J1 \ ( (SecDev(i,:))' - J2 * XYDev(i,:)' );
                        
                        MatrixH1(:, 3*i-2:3*i)  = [dNd2(1) 0 0; 0 dNd2(3) 0; dNd2(3) dNd2(1) 0];
                        MatrixH2(:, 3*i-2:3*i)  = [dNd2(3) 0 0; 0 dNd2(2) 0; dNd2(2) dNd2(3) 0];
                        MatrixG1(:, 3*i-2:3*i)  = [0 0 XYDev(i,1)];
                        MatrixG2(:, 3*i-2:3*i)  = [0 0 XYDev(i,2)];
                        MatrixM1(:, 3*i-2:3*i)  = [XYDev(i,1) 0 0; XYDev(i,2) XYDev(i,1) 0];
                        MatrixM2(:, 3*i-2:3*i)  = [0 0 dNd2(1); 0 0 2*dNd2(3)];
                        MatrixN1(:, 3*i-2:3*i)  = [0 XYDev(i,2) 0];
                        MatrixN2(:, 3*i-2:3*i)  = [0 0 dNd2(2)];
                    end
                    
                    % Calculate PreliminaryStressMatrix in the Current Gauss Point - All Stresses are Linear and Discontinuous at Layer Interfaces.
                    % This is OK for SigmaXX, SigmaYY and TauXY. But Interlaminar Stresses TauXZ TauYZ and SigmaZZ need additional calculation.
                    PreliminaryStressMatrixGaussXYZ([1 2 6 3], 1, layer, gp) =  CBar([1 2 6 3], [1 2 6]) * MatrixB1 * Disp_Bottom  +  CBar([1 2 6 3], 3) * MatrixB2 * (Disp_Top-Disp_Bottom) / h;
                    PreliminaryStressMatrixGaussXYZ([1 2 6 3], 2, layer, gp) =  CBar([1 2 6 3], [1 2 6]) * MatrixB1 * Disp_Top     +  CBar([1 2 6 3], 3) * MatrixB2 * (Disp_Top-Disp_Bottom) / h;
                    PreliminaryStressMatrixGaussXYZ([5 4],     1, layer, gp) =  CBar([5 4], [5 4]) * MatrixBbar2 * Disp_Bottom  +  CBar([5 4], [5 4]) * MatrixBbar1 * (Disp_Top-Disp_Bottom) / h;
                    PreliminaryStressMatrixGaussXYZ([5 4],     2, layer, gp) =  CBar([5 4], [5 4]) * MatrixBbar2 * Disp_Top     +  CBar([5 4], [5 4]) * MatrixBbar1 * (Disp_Top-Disp_Bottom) / h;
                    
                    % Calculate Derivatives dTauXZ/dX dTauYZ/dY dSigmaZZ/dZ in the Current Gauss Point
                    derStressMatrixGaussXYZ([1 2], 1, layer, gp) = -CBar([1 6], [1 2 6]) * MatrixH1 * Disp_Bottom   - CBar([6 2], [1 2 6]) * MatrixH2 * Disp_Bottom   - CBar([1 6], 3) * MatrixG1 * (Disp_Top-Disp_Bottom) / h   -  CBar([6 2], 3) * MatrixG2 * (Disp_Top-Disp_Bottom) / h;
                    derStressMatrixGaussXYZ([1 2], 2, layer, gp) = -CBar([1 6], [1 2 6]) * MatrixH1 * Disp_Top      - CBar([6 2], [1 2 6]) * MatrixH2 * Disp_Top      - CBar([1 6], 3) * MatrixG1 * (Disp_Top-Disp_Bottom) / h   -  CBar([6 2], 3) * MatrixG2 * (Disp_Top-Disp_Bottom) / h;
                    derStressMatrixGaussXYZ(3,     1, layer, gp) = -CBar(5,     [5 4])   * MatrixM2 * Disp_Bottom   - CBar(4, 4)           * MatrixN2 * Disp_Bottom   - CBar(5, [5 4]) * MatrixM1 * (Disp_Top-Disp_Bottom) / h   -  CBar(4, 4)     * MatrixN1 * (Disp_Top-Disp_Bottom) / h;
                    derStressMatrixGaussXYZ(3,     2, layer, gp) = -CBar(5,     [5 4])   * MatrixM2 * Disp_Top      - CBar(4, 4)           * MatrixN2 * Disp_Top      - CBar(5, [5 4]) * MatrixM1 * (Disp_Top-Disp_Bottom) / h   -  CBar(4, 4)     * MatrixN1 * (Disp_Top-Disp_Bottom) / h;
                    
                    TauVectorGaussXYZ(nlayers+layer+1, 1:3, gp) = (PreliminaryStressMatrixGaussXYZ([5 4 3],1, layer, gp) + PreliminaryStressMatrixGaussXYZ([5 4 3],2, layer, gp))/2*h;
                    
                end% Gauss Points
            end%Layers
            
            %Additional Calculation of TauVectorGauss Members
            for layer = 2*nlayers+2 : 3*nlayers
                J = layer - 2*nlayers - 1;
                for gp = 1:ngaus
                    TauVectorGaussXYZ(layer, 1:3, gp) = derStressMatrixGaussXYZ(1:3, 2, J, gp) - derStressMatrixGaussXYZ(1:3, 1, J+1, gp);
                end% Gauss Points
            end%Layers
            
            
            for layer = 1:nlayers
                
                % Extract Thickness of the Current Numerical Layer of the Finite Element
                h = obj.Section.Layers(layer).Thickness;
                
                for gp = 1:ngaus
                    StressMatrixGaussXYZ([1 2 6], [1 2], layer, gp) = PreliminaryStressMatrixGaussXYZ([1 2 6], [1 2], layer, gp);
                    VectorCGaussXYZ(:, 1, gp) = obj.Section.MatrixL \ TauVectorGaussXYZ(:, 1, gp);
                    VectorCGaussXYZ(:, 2, gp) = obj.Section.MatrixL \ TauVectorGaussXYZ(:, 2, gp);
                    VectorCGaussXYZ(:, 3, gp) = obj.Section.MatrixL \ TauVectorGaussXYZ(:, 3, gp);
                    
                    StressMatrixGaussXYZ([5 4], 1, layer, gp) =  VectorCGaussXYZ(3*layer,   [1 2], gp)';
                    StressMatrixGaussXYZ([5 4], 2, layer, gp) =  VectorCGaussXYZ(3*layer-2, [1 2], gp)' * h^2  +  VectorCGaussXYZ(3*layer-1, [1 2], gp)' * h  +  VectorCGaussXYZ(3*layer, [1 2], gp)';
                    
                    StressMatrixGaussXYZ(3, 1, layer, gp) =  VectorCGaussXYZ(3*layer,   3, gp);
                    StressMatrixGaussXYZ(3, 2, layer, gp) =  VectorCGaussXYZ(3*layer-2, 3, gp) * h^2  +  VectorCGaussXYZ(3*layer-1, 3, gp) * h  +  VectorCGaussXYZ(3*layer, 3, gp);
                    
                    
                    % Vector of Constants ai, bi, ci for the Quadratic Interpolation of TauXZ, TauYZ and SigmaZZ
                    % Rows    = ai, bi, ci
                    % Columns = Layers
                    % Third   = TauXZ, TauYZ, SigmaZZ
                    CONSTANTS = zeros(3,nlayers,3);
                    for stress = 1:3
                        CONSTANTS(:, layer, stress) = VectorCGaussXYZ(3*layer-2:3*layer , stress, gp);
                    end
                    
                    for micro = 1:sublayers+1
                        A = StressMatrixGaussXYZ([1 2 6], 1, layer, gp);
                        B = StressMatrixGaussXYZ([1 2 6], 2, layer, gp);
                        obj.StressesXYZ([1 2 6], micro, layer, gp)  =  A + (B-A)/sublayers*(micro-1);
                        
                        SubThickness = (micro-1)*h/sublayers;
                        obj.StressesXYZ(5, micro, layer, gp)  =  CONSTANTS(1,layer,1)*SubThickness^2 + CONSTANTS(2,layer,1)*SubThickness + CONSTANTS(3,layer,1);
                        obj.StressesXYZ(4, micro, layer, gp)  =  CONSTANTS(1,layer,2)*SubThickness^2 + CONSTANTS(2,layer,2)*SubThickness + CONSTANTS(3,layer,2);
                        obj.StressesXYZ(3, micro, layer, gp)  =  CONSTANTS(1,layer,3)*SubThickness^2 + CONSTANTS(2,layer,3)*SubThickness + CONSTANTS(3,layer,3);
                        
                        % Transformation of Stress Components to Material Coordinates of the Current Layer
                        obj.Stresses123(:, micro, layer, gp) = obj.Section.Layers(layer).MatrixR * obj.StressesXYZ(:, micro, layer, gp);
                        
                    end % MicroLayer
                end % Gauss Points
                
            end%Layers
                
        end
      
    end%methods
end