%--------------------------------------------------------------------------
% FLWTFEM MATLAB SOLVER for FINITE ELEMENT ANALYSIS OF MULTILAYER PLATES
% BASED ON THE FULL LAYERWISE THEORY OF REDDY
% Developed by: Miroslav Marjanovic, PhD Civil Eng.
%               Emilija Damnjanovic, MSc Civil Eng.
%               Belgrade, 2019.
%--------------------------------------------------------------------------
function [Points,Weights] = GaussCoordinates(elemshape, Nq)
%--------------------------------------------------------------------------
% This function provides coordinates and weights of 2D Gauss-Legendre
% quadrature points for prescribed element shape and level of numerical
% integration.
%
%  INPUT DATA:
%      elemtype - the shape of 2D finite element (Triangle, Quadrilateral)
%      Nq       - number of quadrature points
%
%  OUTPUT DATA:
%      Points   - coordinates of Gauss points (xi, eta)
%                 size(Points) = Nq x 2
%      Weights  - Weights of Gauss-Legendre quadrature
%                 size(Weights) = Nq x 1
%--------------------------------------------------------------------------
switch elemshape
    
    case 'Quadrilateral'
        switch Nq
            case 1
                Points = [0 0];
                Weights = 4;
                
            case 4
                Points = sqrt(3)/3*[-1 -1; 1  -1; 1 1; -1  1];
                Weights = [1 1 1 1];
                
            case 9
                Points = [-0.774596669241483  -0.774596669241483;...
                           0.774596669241483  -0.774596669241483;...
                           0.774596669241483   0.774596669241483;...
                          -0.774596669241483   0.774596669241483;...
                           0                  -0.774596669241483;...
                           0.774596669241483   0;...
                           0                   0.774596669241483
                          -0.774596669241483   0;...
                           0                   0];
                Weights = 1/81 * [25 25 25 25 40 40 40 40 64]';
        end
        
        
    case 'Triangle'
        switch Nq
                
            case 3
                %Gauss Points are on the Middle of Element Sides
                Points   = [1/2 1/2;
                            0 1/2;
                            1/2 0];
                Weights  =  1/2 * [1/3; 1/3; 1/3];
        end
end