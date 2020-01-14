%--------------------------------------------------------------------------
% FLWTFEM MATLAB SOLVER for FINITE ELEMENT ANALYSIS OF MULTILAYER PLATES
% BASED ON THE FULL LAYERWISE THEORY OF REDDY 
% Developed by: Miroslav Marjanovic, PhD Civil Eng.
%               Emilija Damnjanovic, MSc Civil Eng.
%               Belgrade, 2019.
%--------------------------------------------------------------------------
classdef MathematicalModel < handle
    
    properties
        % Characteristic System Matrices and Vectors
        SystemMatrixKL;
        SystemMatrixM;
        SystemVectorQ;
        SystemVectorP;
        
        % Boundary Conditions
        ZeroDOFs;
        UnknownDOFs;
    end

    methods
        % Calculation of Global System Matrices - LeftHand Side
        function obj = calcSystemMatrices(obj, elements, SDOF)
            elemdofs = length(elements(1).IndexDOF);
            elems    = length(elements);
            
            IG = zeros( elemdofs^2, elems );
            JG = zeros( elemdofs^2, elems );
            KG = zeros( elemdofs^2, elems );
            MG = zeros( elemdofs^2, elems );
            
            for elem = 1:elems
                IG(:,elem) = repmat(elements(elem).IndexDOF,  elemdofs, 1);
                JJ = repmat(elements(elem).IndexDOF', elemdofs, 1);
                JG(:,elem) = JJ(:);
                KG(:,elem) = elements(elem).MatrixKL(:);
                MG(:,elem) = elements(elem).MatrixM(:);
            end
            obj.SystemMatrixKL = sparse(IG(:), JG(:), KG(:), SDOF, SDOF);
            obj.SystemMatrixM  = sparse(IG(:), JG(:), MG(:), SDOF, SDOF);
        end
        
        
        % Calculation of Global System Vectors - RightHand Side
        function obj = calcSystemVectors(obj, nodes, elements, NDOF, SDOF)
            QQ = zeros(SDOF,1);
            PP = zeros(SDOF,1);
            
            for i = 1:length(elements)
                QQ(elements(i).IndexDOF, 1)  =  QQ(elements(i).IndexDOF, 1) + elements(i).VectorQ;
            end
            
            for i = 1:length(nodes)
                PP((i-1)*NDOF+1 : (i-1)*NDOF+3, 1) = [nodes(i).ForceX_bottom  nodes(i).ForceY_bottom  nodes(i).ForceZ_bottom]';
                PP(i*NDOF-2     : i*NDOF      , 1) = [nodes(i).ForceX_top     nodes(i).ForceY_top     nodes(i).ForceZ_top]';
            end
            
            obj.SystemVectorP = sparse(PP);
            obj.SystemVectorQ = sparse(QQ);
        end
        
        
        % Calculation of Zero (Known) and Unknown DOFs
        function obj = deriveZeroDOFs(obj, DOFS)
            knownDOFs = zeros(length(DOFS),1);
            for i = 1:length(DOFS)
                if DOFS(i).Free == 0
                    knownDOFs(i) = i;
                end
            end
            knownDOFs = knownDOFs(knownDOFs~=0);
            
            obj.ZeroDOFs = knownDOFs;
            obj.UnknownDOFs = setdiff( 1:length(DOFS), knownDOFs  );
        end
    end
end