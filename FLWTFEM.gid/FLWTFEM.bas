% ----------------------------------------------------------------------------------
% ANALYSIS DATA

ANALYSISTYPE = '*GenData(Analysis_Type:)';
if strcmp(ANALYSISTYPE,'FreeVibration') == 1
  MODES = *GenData(Number_of_Modes:);
else
  SUBLAYERS = *GenData(SubLayers:);
end



% ----------------------------------------------------------------------------------
% COORDINATES AND CONNECTIVITY DATA

COORDINATES = [
*loop nodes
*format "%12.5f %12.5f %12.5f"
*NodesCoord(1) *NodesCoord(2) *NodesCoord(3); 
*end nodes
];


NELEM = *nelem;
NNODE = *npoin;



% ----------------------------------------------------------------------------------
% GENERATION OF MATERIALS AND LAMINAS DATA

count1 = 0;
count2 = 0;
*loop materials *notused
*if(strcmp(matprop(1),"Orthotropic") == 0)
  count1 = count1 + 1;
  MATERIALS(count1) = OrthotropicMaterial('*matprop(0)', *matprop(2),*matprop(3),*matprop(4), *matprop(5),*matprop(6),*matprop(7), *matprop(8),*matprop(9),*matprop(10), *matprop(11));
*endif
*end materials

*loop materials *notused
*if(strcmp(matprop(1),"Lamina") == 0)
  count2 = count2 + 1;
  for i = 1:count1
    if strcmp('*matprop(2)', MATERIALS(i).MaterialName) == 1
      LAMINAS(count2) = OrthotropicLamina('*matprop(0)', MATERIALS(i),*matprop(3),*matprop(4));
    end
  end    
*endif
*end materials
clear count1 count2



% ----------------------------------------------------------------------------------
% ASSIGNED DATA TO ELEMENTS

ELEMENTSDATA = cell(*nelem,5);
*loop elems
  ELEMENTSDATA{*ElemsNum,5} = [0 0 0 0 0 0];
*end elems

*set cond Composite_Section *elems 
*loop elems    
  Layers = strsplit('*cond(2)');
  LAM = [];
  for i = 1:length(Layers)  
    for lam = 1:length(LAMINAS)
      if strcmp(Layers(i), LAMINAS(lam).LaminaName) == 1
        LAM = horzcat(LAM, LAMINAS(lam));
      end
    end
  end

  ELEMENTSDATA{*ElemsNum,2} = '*cond(3)';
  ELEMENTSDATA{*ElemsNum,3} = CompositeSection('*cond(1)', LAM); 
*end elems
clear Layers LAM lam i

*loop elems
*if(ElemsType==2)
*if(ElemsNnode==3)
*format "%6i %6i %6i"
  ELEMENTSDATA{*ElemsNum,4} = [*ElemsConec];
  ELEMENTSDATA{*ElemsNum,1} = 'T3';
*endif
*if(ElemsNnode==6)
*format "%6i %6i %6i %6i %6i %6i"
  ELEMENTSDATA{*ElemsNum,4} = [*ElemsConec];
  ELEMENTSDATA{*ElemsNum,1} = 'T6';
*endif
*endif
*if(ElemsType==3)
*if(ElemsNnode==4)
*format "%6i %6i %6i %6i"
  ELEMENTSDATA{*ElemsNum,4} = [*ElemsConec];
  ELEMENTSDATA{*ElemsNum,1} = 'Q4';
*endif
*if(ElemsNnode==8)
*format "%6i %6i %6i %6i %6i %6i %6i %6i"		
  ELEMENTSDATA{*ElemsNum,4} = [*ElemsConec];
  ELEMENTSDATA{*ElemsNum,1} = 'Q8';
*endif
*endif
*end elems

*set cond Distributed_Loadings *elems 
*loop elems *OnlyInCond
*format "%6i %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f"
  ELEMENTSDATA{*ElemsNum,5} = [*cond(1) *cond(2) *cond(3) *cond(4) *cond(5) *cond(6)];
*end elems



% ----------------------------------------------------------------------------------
% ASSIGNED DATA TO POINTS

POINTSDATA = cell(*npoin,2);
*loop nodes
  POINTSDATA{*NodesNum,1} = [0 0 0 0 0 0];
  POINTSDATA{*NodesNum,2} = [0 0 0 0 0 0 0 0 0];
*end nodes

*set cond Nodal_Forces *nodes
*loop nodes *OnlyInCond
*format "%6i %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f"
  POINTSDATA{*NodesNum,1} = [*cond(1) *cond(2) *cond(3) *cond(4) *cond(5) *cond(6)];
*end nodes

*set cond Point_Constraints *nodes
*loop nodes *OnlyInCond
*format "%6i %6i %6i %6i %6i %6i %6i %6i %6i %6i"
  POINTSDATA{*NodesNum,2} = [*cond(1) *cond(2) *cond(3) *cond(4) *cond(5) *cond(6) *cond(7) *cond(8) *cond(9)];
*end nodes


% ----------------------------------------------------------------------------------
% ADDITIONAL CALCULATIONS and CLEARING
LAYERS = length(ELEMENTSDATA{1,3}.Layers);
NDOF = 3 ** (LAYERS + 1);
SDOF = NNODE ** NDOF;


% ----------------------------------------------------------------------------------
% CALCULATION OF SECTION STIFFNESS AND INERTIA COEFFICIENTS
for elem = 1:size(ELEMENTSDATA,1)
    ELEMENTSDATA{elem,3} = ELEMENTSDATA{elem,3}.calcCoeffMatrices();
    ELEMENTSDATA{elem,3} = ELEMENTSDATA{elem,3}.calcInertiaMatrices();
end
clear elem
clear LAMINAS MATERIALS