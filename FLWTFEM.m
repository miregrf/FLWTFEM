profile on
clc; clear;
%--------------------------------------------------------------------------
% FLWTFEM MATLAB SOLVER for FINITE ELEMENT ANALYSIS OF MULTILAYER PLATES
% BASED ON THE FULL LAYERWISE THEORY OF REDDY 
% Developed by: Miroslav Marjanovic, PhD Civil Eng.
%               Emilija Damnjanovic, MSc Civil Eng.
%               Belgrade, 2019.
%--------------------------------------------------------------------------
%%  EVALUATION of INPUTFILE
disp('READING INPUTFILE...')

    % Add Paths   
    addpath([pwd '\materialClasses\']);  addpath([pwd '\modelClasses\']);
    addpath([pwd '\utilityFunctions\']); addpath([pwd '\postAlgorithms\']);
    
    % Read GiD INPUT FILE
    DIR = uigetdir;                      [~,NAME] = fileparts(DIR);
    source = [DIR, '\', NAME,'.dat'];    copyfile(source,pwd,'f');
    FILENAME = [NAME,'.dat'];            INPUTFILE = fileread(FILENAME);
    eval(INPUTFILE);                     
    clear source INPUTFILE
    tic
%--------------------------------------------------------------------------
%  After reading of the INPUTFILE, the following variables are created:

%  ANALYSISTYPE = Type of Analysis (Bending or FreeVibration)
%  NELEM = Total Number of Finite Elements
%  NNODE = Total Number of Nodes
%  COORDINATES  = NNODE x 3 Matrix of Nodal Coordinates (X, Y, Z)

%  ELEMENTSDATA = NELEM x 5 Cell with the Following Data in Columns:
%    1 - ELEMENT TYPE
%    2 - TYPE OF INTEGRATION
%    3 - ASSIGNED SECTION
%    4 - CONNECTIVITY VECTOR
%    5 - ASSIGNED DISTRIBUTED LOADINGS
%  POINTSDATA = NNODE x 2 CELL WITH THE FOLLOWING DATA IN COLUMNS:
%    1 - ASSIGNED NODAL FORCES
%    2 - ASSIGNED POINT CONSTRAINTS

%  LAYERS = Number of Material Layers in the Considered Plate
%  NDOF   = Number of Degrees Of Freedom in the Node
%  SDOF   = Total Number of Degrees Of Freedom in the Model

%  MODES = Number of Calculated EigenModes in FreeVibration Analysis
%  SUBLAYERS = Number of SubLayers for Interlaminar Stress Calculation in Bending Analysis

%  DIR      = Model Directory
%  FILENAME = Model File Name
%  NAME     = Model Name
%--------------------------------------------------------------------------    
%% GENERATION OF GLOBAL COMPUTATIONAL MODEL
disp('GENERATION OF GLOBAL COMPUTATIONAL MODEL...')
    
    % Construction of the MODEL Instance of MathematicalModel Class
    MODEL = MathematicalModel();
    
    % Construction of DOFS Instances of DOF Class
    DOFS(SDOF) = DOF();
    for dof = 1:SDOF
        DOFS(dof) = DOF(dof);
    end
    clear dof
    
    
%%  GENERATION OF NODES, ASSIGNMENT OF DOFs, ASSIGNMENT OF NODAL FORCES
disp('GENERATION OF NODES, ASSIGNMENT OF DOFs AND NODAL FORCES...')

    NODES(NNODE) = Node();
    for node = 1:NNODE
        % Construction of NODES Instances of Node Class
        NODES(node) = Node(node,COORDINATES);
        % Assignment of DOFs and Nodal Forces to Nodes
        NODES(node) = NODES(node).assignDOFs(NDOF, DOFS);
        NODES(node) = NODES(node).assignNodalForces(POINTSDATA{node,1});
    end
    clear node COORDINATES
    
    
%% ASIGNMENT OF BOUNDARY CONDITIONS
disp('ASSIGNMENT OF BOUNDARY CONDITIONS...')

    % Assignment of Boundary Conditions to DOFS and NODES
    DOFS = ConstraintHandler(DOFS, NODES, POINTSDATA);
    for node = 1:NNODE
        NODES(node) = NODES(node).assignDOFs(NDOF, DOFS);
    end
    
    MODEL = MODEL.deriveZeroDOFs(DOFS);
    clear node POINTSDATA DOFS


%%  GENERATION OF FINITE ELEMENTS and CALCULATION of CHARACTERISTIC MATRICES
disp('GENERATION OF FINITE ELEMENTS...')
    
    FEs = LW_3D.empty(NELEM,0);
    for elem = 1:NELEM
        % Construction of FEs Instances of LW_3D<LW_L1 Class
        FEs(elem) = LW_3D(elem, ELEMENTSDATA{elem,1}, ELEMENTSDATA{elem,2}, ELEMENTSDATA{elem,3}, NODES(ELEMENTSDATA{elem,4}), ELEMENTSDATA{elem,5});
        
        % Calculation of Characteristic Element Matrices and Vectors
        FEs(elem) = FEs(elem).deriveGauss();
        FEs(elem) = FEs(elem).calcMatrixKL();
        FEs(elem) = FEs(elem).calcMatrixM();
        FEs(elem) = FEs(elem).calcVectorQ();
        FEs(elem) = FEs(elem).deriveIndexDOF();
    end
    clear elem node ELEMENTSDATA


%% GENERATION OF GLOBAL COMPUTATIONAL MODEL
disp('GENERATION OF GLOBAL COMPUTATIONAL MODEL...')
    MODEL = MODEL.calcSystemMatrices(FEs, SDOF);
    MODEL = MODEL.calcSystemVectors(NODES, FEs, NDOF, SDOF);
   
   
%% ANALYSIS
    Knn = MODEL.SystemMatrixKL(MODEL.UnknownDOFs, MODEL.UnknownDOFs);

    % FREE VIBRATION ANALYSIS
    if strcmp(ANALYSISTYPE,'FreeVibration') == 1
        SUBLAYERS = 1;
        Mnn = MODEL.SystemMatrixM(MODEL.UnknownDOFs, MODEL.UnknownDOFs);
        
        [vector, values] = eigs(Knn, Mnn, MODES, 0);
        values = diag(sqrt(values));
        [values,number] = sort(values);
        number = number(1 : MODES);
        NatFreqs = values(1:MODES);

        FINALDISP = zeros(SDOF, MODES);
        for mode = 1:MODES
            FINALDISP(MODEL.UnknownDOFs, mode) = vector(:,number(mode));
        end

        % Assignment of Results to Element Instances
        for elem = 1:NELEM
            FEs(elem).Displacements = FINALDISP( FEs(elem).IndexDOF,  :);
        end
    end

    % BENDING ANALYSIS
    if strcmp(ANALYSISTYPE,'Bending') == 1
        Sn  = MODEL.SystemVectorQ(MODEL.UnknownDOFs, 1) + MODEL.SystemVectorP(MODEL.UnknownDOFs, 1);

        FINALDISP = zeros(SDOF, 1);
        FINALDISP(MODEL.UnknownDOFs, 1)  =  Knn \ Sn;

        for elem = 1:NELEM
            % Assignment of Results to Element Instances 
            FEs(elem).Displacements = FINALDISP( FEs(elem).IndexDOF,  1);
            % Post-Processing of Stresses
            FEs(elem) = FEs(elem).CalcStresses(SUBLAYERS);
        end
    end
    clear Knn Mnn Sn elem vector values number mode
    

%% GENERATION OF POST PROCESS NODES
disp('GENERATION OF POST PROCESS NODES and ELEMENTS...')
   
    % Calculation of z-coordinate of Every SubLayer along Plate Thickness
    zz = zeros(LAYERS*SUBLAYERS+1, 1);
    hh = 0;
    count = 1;
    for layer = 1:LAYERS
        for sublayer = 1:SUBLAYERS
            hh = hh + FEs(1).Section.Layers(layer).Thickness / SUBLAYERS;
            count = count + 1;
            zz(count) = hh;
        end
    end
    zz = zz - hh/2;
    
    % Node Indices for the New PostNode Instances are first assigned in the
    % first (bottom) numerical layer. After that, starting from the "first"
    % node in the 2D mesh, the indices are assigned in the next layer going
    % from the bottom to the top of the plate.
    POSTNODES( NNODE * (LAYERS*SUBLAYERS + 1) ) = PostNode();
    for layer = 1:LAYERS*SUBLAYERS + 1
        for node = 1:NNODE
            % Construction of POSTNODES Instances of POSTNODE Class
            POSTNODES(  (layer-1)*NNODE + node  ) = PostNode((layer-1)*NNODE + node,    NODES(node).XCoord, NODES(node).YCoord, zz(layer));
        end
    end
    clear zz hh sublayer count layer node
    
    for node = 1:NNODE
        % Calculated Displacements in the Current Node of the 2D Mesh
        CurrentNodeDisp = FINALDISP( (node-1)*NDOF+1 : node*NDOF,  :);
        
        % Assignment of Results in the first layer of Post-Process Nodes
        POSTNODES( node ).DispX = CurrentNodeDisp(1, :);
        POSTNODES( node ).DispY = CurrentNodeDisp(2, :);
        POSTNODES( node ).DispZ = CurrentNodeDisp(3, :);
                
        count = 0;
        for layer = 1:LAYERS
            X_bottom = CurrentNodeDisp(3*layer-2, :);
            X_top    = CurrentNodeDisp(3*layer+1, :);
            Y_bottom = CurrentNodeDisp(3*layer-1, :);
            Y_top    = CurrentNodeDisp(3*layer+2, :);
            Z_bottom = CurrentNodeDisp(3*layer  , :);
            Z_top    = CurrentNodeDisp(3*layer+3, :);
            
            for sublayer = 1:SUBLAYERS
                count = count + 1;
                % Assignment of Results in all layers of Post-Process
                % Nodes. Linear Interpolation of Displacements is done.
                POSTNODES(  count*NNODE + node  ).DispX = X_bottom + (X_top-X_bottom) / SUBLAYERS * sublayer;
                POSTNODES(  count*NNODE + node  ).DispY = Y_bottom + (Y_top-Y_bottom) / SUBLAYERS * sublayer;
                POSTNODES(  count*NNODE + node  ).DispZ = Z_bottom + (Z_top-Z_bottom) / SUBLAYERS * sublayer;
            end
        
        end
    end
    clear node CurrentNodeDisp count layer sublayer
    clear X_bottom Y_bottom Z_bottom X_top Y_top Z_top
    

%% GENERATION OF POST PROCESS ELEMENTS    
    count = 0;
    countlayer = 0;
    POSTELEMENTS(NELEM*LAYERS*SUBLAYERS) = PostElement();
    for layer = 1:LAYERS
        for sublayer = 1:SUBLAYERS
            countlayer = countlayer+1;
            for elem = 1:NELEM
                
                if strcmp(FEs(elem).ElementType, 'T3') == 1 || strcmp(FEs(elem).ElementType, 'T6') == 1
                    % For Both T3 and T6, the GiD Output will be 6-Node Prism, Including only the Corner Nodes of the 2D Element
                    array = (countlayer-1)*NNODE + [FEs(elem).ElementNodes(1:3).NodeID];
                elseif strcmp(FEs(elem).ElementType, 'Q4') == 1 || strcmp(FEs(elem).ElementType, 'Q8') == 1
                    % For Both Q4 and Q8, the GiD Output will be 8-Node Prism, Including only the Corner Nodes of the 2D Element
                    array = (countlayer-1)*NNODE + [FEs(elem).ElementNodes(1:4).NodeID];
                end
                
                count = count + 1;
                % Construction of POSTELEMENTS Instances of POSTELEMENT Class
                POSTELEMENTS(count) = PostElement(count, POSTNODES([array array+NNODE]));
                
                
                if strcmp(ANALYSISTYPE,'Bending') == 1

                    ngaus    = size(FEs(elem).StressesXYZ, 4);
                    
                    % POSTELEMENTS.StressesXYZ and POSTELEMENTS.Stresses123
                    % Rows    = Six Componental Stresses SigmaXX SigmaYY SigmaZZ TauYZ TauXZ TauXY
                    % Columns = Two Interfaces (Bottom and Top) of the Current Elem (Layer)
                    % Third   = Gauss Points
                    POSTELEMENTS(count).StressesXYZ = zeros(6,2,ngaus);
                    POSTELEMENTS(count).Stresses123 = zeros(6,2,ngaus);
                    
                    for gaus = 1:ngaus
                        POSTELEMENTS(count).StressesXYZ(:,:,gaus) = FEs(elem).StressesXYZ(:,[sublayer sublayer+1],layer,gaus);
                        POSTELEMENTS(count).Stresses123(:,:,gaus) = FEs(elem).Stresses123(:,[sublayer sublayer+1],layer,gaus);
                    end
                    
                end
            end
        end
    end
    clear count countlayer layer sublayer elem array
    clear gaus ngaus


%%  GID OUTPUT
disp('GENERATING GiD OUTPUT...')
    make_post_msh(FILENAME, POSTELEMENTS, POSTNODES);

    if strcmp(ANALYSISTYPE,'Bending') == 1   
        make_post_res(FILENAME, POSTNODES, POSTELEMENTS, 'Bending', 1)    
    else
        make_post_res(FILENAME, POSTNODES, POSTELEMENTS, 'FreeVibration', NatFreqs)    
    end
    

%%  CLOSE ALL
	close all
	movefile([pwd,'\', NAME, '.post.msh'],  DIR,  'f');
    movefile([pwd,'\', NAME, '.post.res'],  DIR,  'f');
    movefile([pwd,'\', NAME, '.dat'],       DIR,  'f');

    rmpath([pwd '\materialClasses\']);   rmpath([pwd '\modelClasses\']);
    rmpath([pwd '\utilityFunctions\']);  rmpath([pwd '\postAlgorithms\']);
    
    clear name DIR FILENAME;    toc