%--------------------------------------------------------------------------
% FLWTFEM MATLAB SOLVER for FINITE ELEMENT ANALYSIS OF MULTILAYER PLATES
% BASED ON THE FULL LAYERWISE THEORY OF REDDY 
% Developed by: Miroslav Marjanovic, PhD Civil Eng.
%               Emilija Damnjanovic, MSc Civil Eng.
%               Belgrade, 2019.
%--------------------------------------------------------------------------
function DOFS = ConstraintHandler(DOFS, NODES, POINTSDATA)
    NNODE = length(NODES);

    for node = 1:NNODE
        % Assign Particular Nodal Restraints to DOFs
        bc = POINTSDATA{node,2};
        %Ux_top
        if bc(1) == 1
            DOFS(NODES(node).NodalDOFs(end-2).ID).Free = 0;
        end
        %Uy_top
        if bc(2) == 1
            DOFS(NODES(node).NodalDOFs(end-1).ID).Free = 0;
        end
        %Uz_top
        if bc(3) == 1
            DOFS(NODES(node).NodalDOFs(end).ID).Free = 0;
        end
        %Ux_bottom
        if bc(4) == 1
            DOFS(NODES(node).NodalDOFs(1).ID).Free = 0;
        end
        %Uy_bottom
        if bc(5) == 1
            DOFS(NODES(node).NodalDOFs(2).ID).Free = 0;
        end
        %Uz_bottom
        if bc(6) == 1
            DOFS(NODES(node).NodalDOFs(3).ID).Free = 0;
        end
    end

    for node = 1:NNODE
        % Assign Nodal Restraints to all DOFs along the Thickness
        bc = POINTSDATA{node,2};
        %Ux_ALL
        if bc(7) == 1
            array = [NODES(node).NodalDOFs(1:3:end-2).ID];
            for i = 1:length(array)
                DOFS(array(i)).Free = 0;
            end
        end
        %Uy_ALL
        if bc(8) == 1
            array = [NODES(node).NodalDOFs(2:3:end-1).ID];
            for i = 1:length(array)
                DOFS(array(i)).Free = 0;
            end
        end
        %Uz_ALL
        if bc(9) == 1
            array = [NODES(node).NodalDOFs(3:3:end).ID];
            for i = 1:length(array)
                DOFS(array(i)).Free = 0;
            end
        end
    end