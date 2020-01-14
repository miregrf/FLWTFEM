# ----------------------------------------------------------------------------
# LayerWiseFE.tcl  *.tcl Script program for the LayerWiseFE GiD problemtype
# Ruhr University Bochum
# Miroslav Marjanovic, M.Sc. B.Sc. Civil Eng.
# ----------------------------------------------------------------------------
# This file is written in TCL language 
# For more information about TCL look at: http://www.sunlabs.com/research/tcl/
# For more information about GID internals, check the program scripts.
#
# At least two procs must be in this file:
# InitGIDProject dir - Will be called whenever a project is begun to be used,
#                      where dir is the project's directory
#
# EndGIDProject - Will be called whenever a project ends to be used.
# ----------------------------------------------------------------------------
#
# prefix values:
#     Pre        Only active in the preprocessor
#     Post       Only active in the postprocessor
#     PrePost    Active Always
#
# ----------------------------------------------------------------------------

proc MyBitmaps { dir { type "DEFAULT INSIDELEFT"} } {
}

proc InitGIDProject { dir } {
    
    GiD_DataBehaviour materials "Orthotropic_Lamina" hide { assign draw unassign delete impexp}
    
    GiDMenu::RemoveOption "Geometry" [list "Create" "Volume"] PRE _
    GiDMenu::RemoveOption "Geometry" [list "Create" "Contact"] PRE _
    GiDMenu::RemoveOption "Geometry" [list "Create" "Object"] PRE _
    GiDMenu::RemoveOption "Geometry" [list "Create" "---3"] PRE _
    GiDMenu::RemoveOption "Geometry" [list "Create" "---2"] PRE _
    
    GiDMenu::RemoveOption "Utilities" [list "Signal"] PREPOST _
    GiDMenu::RemoveOption "Utilities" [list "Swap normals"] PREPOST _
    GiDMenu::RemoveOption "Utilities" [list "Renumber"] PRE _
    GiDMenu::RemoveOption "Utilities" [list "Id"] PREPOST _

    GiDMenu::RemoveOption "Data" [list "Interval"] PRE _
    
    GiDMenu::RemoveOption "Mesh" [list "Unstructured" "Assign sizes on points"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Unstructured" "Assign sizes on volumes"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Unstructured" "Assign entities" "Volumes"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Unstructured" "Sizes by chordal error..."] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Unstructured" "Sizes by background mesh..."] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Unstructured" "---0"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Unstructured" "---0"] PRE 
    GiDMenu::RemoveOption "Mesh" [list "Structured" "Volumes"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "SemiStructured"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "---0"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "---0"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "---1"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "---2"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Cartesian"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Boundary layer"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Element type" "Default"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Element type" "Linear"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Element type" "---0"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Element type" "---0"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Element type" "Circle"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Element type" "Tetrahedra"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Element type" "Hexahedra"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Element type" "Prism"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Element type" "Sphere"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Element type" "Only points"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Draw"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Draw" "Skip entities (Rjump)"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Draw" "Force points to"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Draw" "Boundary layer"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Mesh criteria"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Reset mesh data"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Edit mesh"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Create boundary mesh"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Mesh options from model"] PRE _
    GiDMenu::RemoveOption "Mesh" [list "Mesh quality..."] PRE _
    GiDMenu::RemoveOption "Mesh" [list "---1"] PRE _
    
    GiDMenu::RemoveOption "Calculate" [list "Calculate remote"] PRE _
    GiDMenu::RemoveOption "Calculate" [list "View process info..."] PRE _
    GiDMenu::RemoveOption "Calculate" [list "Cancel process"] PRE _
    GiDMenu::Delete "Help" PREPOST _
    
    GiDMenu::UpdateMenus  
}

proc EndGIDProject {} {
}