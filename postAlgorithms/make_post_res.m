%--------------------------------------------------------------------------
% FLWTFEM MATLAB SOLVER for FINITE ELEMENT ANALYSIS OF MULTILAYER PLATES
% BASED ON THE FULL LAYERWISE THEORY OF REDDY 
% Developed by: Miroslav Marjanovic, PhD Civil Eng.
%               Emilija Damnjanovic, MSc Civil Eng.
%               Belgrade, 2019.
%--------------------------------------------------------------------------
function make_post_res(filename, POSTNODES, POSTELEMENTS, ANALYSISTYPE, RECORDS)
    
    %INITIALIZATION
    res_file = strcat(filename(1:end-4),'.post.res');
    fid = fopen(res_file,'w');
    fprintf(fid,'GiD Post Results File 1.0                                   \n');
    
%% GAUSS POINTS COORDINATES
    
    Prism6     = 0;
    Hexahedra8 = 0;
    
    for elem = 1:length(POSTELEMENTS)
        if strcmp(POSTELEMENTS(elem).ElementType, 'Prism 6') == 1
            Prism6 = 1;
            
            fprintf(fid,'GaussPoints "Gauss_Points_Prism_6" Elemtype Prism          \n');
            fprintf(fid,'Number of Gauss Points: 6                                  \n');
            fprintf(fid,'Natural Coordinates: Given                                 \n');
            fprintf(fid,'0.50 0.50 0.00                                             \n');
            fprintf(fid,'0.00 0.50 0.00                                             \n');
            fprintf(fid,'0.50 0.00 0.00                                             \n');
            fprintf(fid,'0.50 0.50 1.00                                             \n');
            fprintf(fid,'0.00 0.50 1.00                                             \n');
            fprintf(fid,'0.50 0.00 1.00                                             \n');
            fprintf(fid,'end GaussPoints                                            \n');
            
            break;
            
        elseif strcmp(POSTELEMENTS(elem).ElementType, 'Hexahedra 8') == 1
            Hexahedra8 = 1;
            
            fprintf(fid,'GaussPoints "Gauss_Points_Hexahedra_8" Elemtype Hexahedra  \n');
            fprintf(fid,'Number of Gauss Points: 8                                  \n');
            fprintf(fid,'Natural Coordinates: Given                                 \n');
            fprintf(fid,'-0.577350269189626 -0.577350269189626 -1                   \n');
            fprintf(fid,' 0.577350269189626 -0.577350269189626 -1                   \n');
            fprintf(fid,' 0.577350269189626  0.577350269189626 -1                   \n');
            fprintf(fid,'-0.577350269189626  0.577350269189626 -1                   \n');
            fprintf(fid,'-0.577350269189626 -0.577350269189626  1                   \n');
            fprintf(fid,' 0.577350269189626 -0.577350269189626  1                   \n');
            fprintf(fid,' 0.577350269189626  0.577350269189626  1                   \n');
            fprintf(fid,'-0.577350269189626  0.577350269189626  1                   \n');
            fprintf(fid,'end GaussPoints                                            \n');
            
            break;

        end
    end
    
    
%% BENDING ANALYSIS
    if strcmp(ANALYSISTYPE, 'Bending') == 1

        % Deformed Shape
        fprintf(fid,'Result "Deformed Shape" "Bending Analysis" %6i Vector OnNodes \n', 1);
        fprintf(fid,'ComponentNames "X", "Y", "Z" \n');

        fprintf(fid,'Values \n');

        fprintf(fid,'%6i %12.5f %12.5f %12.5f \n', [[1:size(POSTNODES, 2)]; POSTNODES.DispX; POSTNODES.DispY; POSTNODES.DispZ]);

        fprintf(fid,'End Values \n');

        
        % Stresses
        % Prism 6
        if Prism6 == 1
            
            % Stresses XYZ
            fprintf(fid,'Result "Stresses XYZ" "Bending Analysis" %6i Matrix OnGaussPoints "Gauss_Points_Prism_6"  \n', 1);
            fprintf(fid,'ComponentNames "SigmaXX", "SigmaYY", "SigmaZZ", "SigmaYZ", "SigmaXZ", "SigmaXY" \n');
            fprintf(fid,'Values \n');
            for elem = 1:length(POSTELEMENTS)
                if strcmp(POSTELEMENTS(elem).ElementType, 'Prism 6') == 1
                    fprintf(fid,'%6i %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f \n', elem, POSTELEMENTS(elem).StressesXYZ(:, 1, 1)');
                    for i = 2:3
                        fprintf(fid,'    %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f \n',   POSTELEMENTS(elem).StressesXYZ(:, 1, i)');
                    end
                    for i = 1:3
                        fprintf(fid,'    %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f \n',   POSTELEMENTS(elem).StressesXYZ(:, 2, i)');
                    end
                end
            end
            fprintf(fid,'End Values \n');

            fprintf(fid,'Result "Stresses 123" "Bending Analysis" %6i Matrix OnGaussPoints "Gauss_Points_Prism_6"  \n', 1);
            fprintf(fid,'ComponentNames "Sigma1", "Sigma2", "Sigma3", "Sigma4", "Sigma5", "Sigma6" \n');
            fprintf(fid,'Values \n');
            for elem = 1:length(POSTELEMENTS)
                if strcmp(POSTELEMENTS(elem).ElementType, 'Prism 6') == 1
                    fprintf(fid,'%6i %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f \n', elem, POSTELEMENTS(elem).Stresses123(:, 1, 1)');
                    for i = 2:3
                        fprintf(fid,'    %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f \n',   POSTELEMENTS(elem).Stresses123(:, 1, i)');
                    end
                    for i = 1:3
                        fprintf(fid,'    %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f \n',   POSTELEMENTS(elem).Stresses123(:, 2, i)');
                    end
                end
            end
            fprintf(fid,'End Values \n');
        end

        %Hexahedra 8
        if Hexahedra8 == 1
            
            % Stresses XYZ
            fprintf(fid,'Result "Stresses XYZ" "Bending Analysis" %6i Matrix OnGaussPoints "Gauss_Points_Hexahedra_8"  \n', 1);
            fprintf(fid,'ComponentNames "SigmaXX", "SigmaYY", "SigmaZZ", "SigmaYZ", "SigmaXZ", "SigmaXY" \n');
            fprintf(fid,'Values \n');
            for elem = 1:length(POSTELEMENTS)
                if strcmp(POSTELEMENTS(elem).ElementType, 'Hexahedra 8') == 1
                    fprintf(fid,'%6i %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f \n', elem, POSTELEMENTS(elem).StressesXYZ(:, 1, 1)');
                    for i = 2:4
                        fprintf(fid,'    %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f \n',   POSTELEMENTS(elem).StressesXYZ(:, 1, i)');
                    end
                    for i = 1:4
                        fprintf(fid,'    %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f \n',   POSTELEMENTS(elem).StressesXYZ(:, 2, i)');
                    end
                end
            end
            fprintf(fid,'End Values \n');

            % Stresses 123
            fprintf(fid,'Result "Stresses 123" "Bending Analysis" %6i Matrix OnGaussPoints "Gauss_Points_Hexahedra_8"  \n', 1);
            fprintf(fid,'ComponentNames "Sigma1", "Sigma2", "Sigma3", "Sigma4", "Sigma5", "Sigma6" \n');
            fprintf(fid,'Values \n');
            for elem = 1:length(POSTELEMENTS)
                if strcmp(POSTELEMENTS(elem).ElementType, 'Hexahedra 8') == 1
                    fprintf(fid,'%6i %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f \n', elem, POSTELEMENTS(elem).Stresses123([1 2 3 6 5 4], 1, 1)');
                    for i = 2:4
                        fprintf(fid,'    %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f \n',   POSTELEMENTS(elem).Stresses123([1 2 3 6 5 4], 1, i)');
                    end
                    for i = 1:4
                        fprintf(fid,'    %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f \n',   POSTELEMENTS(elem).Stresses123([1 2 3 6 5 4], 2, i)');
                    end
                end
            end
            fprintf(fid,'End Values \n');

        end
    end
        
    
%% FREE VIBRATIONS
    if strcmp(ANALYSISTYPE, 'FreeVibration') == 1
        for record = 1:length(RECORDS)
            
            % Deformed (Mode) Shapes
            fprintf(fid,'Result "Mode Shape" "Eigen Value Analysis - Hz" %6i Vector OnNodes \n', RECORDS(record)/2/pi);
            fprintf(fid,'ComponentNames "X", "Y", "Z" \n');

            fprintf(fid,'Values \n');
            for node = 1:length(POSTNODES)
                fprintf(fid,'%6i %12.5f %12.5f %12.5f \n', node, POSTNODES(node).DispX(record), POSTNODES(node).DispY(record), POSTNODES(node).DispZ(record));
            end
            fprintf(fid,'End Values \n');
        end
    end

    fclose(fid);
end