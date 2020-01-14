%--------------------------------------------------------------------------
% FLWTFEM MATLAB SOLVER for FINITE ELEMENT ANALYSIS OF MULTILAYER PLATES
% BASED ON THE FULL LAYERWISE THEORY OF REDDY 
% Developed by: Miroslav Marjanovic, PhD Civil Eng.
%               Emilija Damnjanovic, MSc Civil Eng.
%               Belgrade, 2019.
%--------------------------------------------------------------------------
function make_post_msh(file_name, POSTELEMENTS, POSTNODES)

    msh_file = strcat(file_name(1:end-4),'.post.msh');
    fid = fopen(msh_file,'w');

        % Make New Mesh of POSTNODES
        fprintf(fid,'MESH dimension 3 elemtype Prism nnode 6 \n');
        fprintf(fid,'coordinates \n');
        for node = 1:length(POSTNODES)
            fprintf(fid,'%6i %12.5f %12.5f %12.5f \n', node, POSTNODES(node).XCoord, POSTNODES(node).YCoord, POSTNODES(node).ZCoord);
        end
        fprintf(fid,'end coordinates \n');
        
        %Make New Connectivity for PRISM Elements according to GiD Customization
        fprintf(fid,'elements \n');
        for elem = 1:length(POSTELEMENTS)
            if length(POSTELEMENTS(elem).ElementNodes) == 6
                fprintf(fid,'%6i %6i %6i %6i %6i %6i %6i \n', elem, [POSTELEMENTS(elem).ElementNodes.NodeID]);
            end
        end
        fprintf(fid,'end elements \n');
        
        %Make New Connectivity for HEXAHEDRA Elements according to GiD Customization
        fprintf(fid,'MESH dimension 3 elemtype Hexahedra nnode 8 \n');
        fprintf(fid,'elements \n');
        for elem = 1:length(POSTELEMENTS)
            if length(POSTELEMENTS(elem).ElementNodes) == 8
                fprintf(fid,'%6i %6i %6i %6i %6i %6i %6i %6i %6i \n', elem, [POSTELEMENTS(elem).ElementNodes.NodeID]);
            end
        end
        fprintf(fid,'end elements');
            
    fclose(fid);

end