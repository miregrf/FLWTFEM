%--------------------------------------------------------------------------
% FLWTFEM MATLAB SOLVER for FINITE ELEMENT ANALYSIS OF MULTILAYER PLATES
% BASED ON THE FULL LAYERWISE THEORY OF REDDY 
% Developed by: Miroslav Marjanovic, PhD Civil Eng.
%               Emilija Damnjanovic, MSc Civil Eng.
%               Belgrade, 2019.
%--------------------------------------------------------------------------
function PlotThicknessNodeStresses(FE, node, cs, LAYERS, SUBLAYERS)

    PlotX1 = zeros(LAYERS*(SUBLAYERS+1), 6);
    PlotY1 = zeros(LAYERS*(SUBLAYERS+1), 1);

    zz       = 0;
    h_all    = 0;
    h_before = 0;
    count    = 0;
    for layer = 1:LAYERS
        h     = FE.Section.Layers(layer).Thickness;
        h_all = h_all + h;
        
        for micro = 1:SUBLAYERS+1
            count = count + 1;
            zz    = zz + h/SUBLAYERS;
            
            if strcmp(cs,'XYZ') == 1
                PlotX1(count,:) = FE.StressesNodalXYZ(:,micro,layer,node)';
            else
                PlotX1(count,:) = FE.StressesNodalXYZ(:,micro,layer,node)';
            end
            PlotY1(count)   = h_before + h/SUBLAYERS*(micro-1);
        end
        
        h_before = h_before+h;
        
    end
    PlotY1 = (PlotY1-h_all/2) / h_all;
    PlotX1(:,4)
    PlotY1
    
%% PLOT

    maxX = max ( max(PlotX1(:,:)) );
    minX = min ( min(PlotX1(:,:)) );
    
    figure
%     suptitle('Stress Distribution through the Thickness of the Plate - NODE')
    plot1 = subplot(2,3,1);
    plot(PlotX1(:,1),PlotY1)
    xlabel('SigmaXX'); ylabel('z / h');
    grid on
%     set(gca,'XLim',[minX, maxX])
    
    plot2 = subplot(2,3,2);
    plot(PlotX1(:,2),PlotY1)
    xlabel('SigmaYY'); ylabel('z / h');
    grid on
% 	set(gca,'XLim',[minX, maxX])
    
    plot3 = subplot(2,3,3);
    plot(PlotX1(:,3),PlotY1)
    xlabel('SigmaZZ'); ylabel('z / h');
    grid on
% 	set(gca,'XLim',[minX, maxX])
    

    plot4 = subplot(2,3,4);
    plot(PlotX1(:,4),PlotY1)
    xlabel('TauYZ'); ylabel('z / h');
    grid on
% 	set(gca,'XLim',[minX, maxX])
    
    plot5 = subplot(2,3,5);
    plot(PlotX1(:,5),PlotY1)
    xlabel('TauXZ'); ylabel('z / h');
    grid on
% 	set(gca,'XLim',[minX, maxX])
    
    plot6 = subplot(2,3,6);
    plot(PlotX1(:,6),PlotY1)
    xlabel('TauXY'); ylabel('z / h');
    grid on
% 	set(gca,'XLim',[minX, maxX])
    
%     linkaxes([plot1,plot2,plot3,plot4,plot5,plot6],'xy')
    
end