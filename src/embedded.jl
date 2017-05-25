
##############################################################################
#    FemLab - Finite Element Library                                         #
#    Copyright (C) 2014 Raul Durand <raul.durand at gmail.com>               #
#                                                                            #
#    This file is part of FemLab.                                            #
#                                                                            #
#    FemLab is free software: you can redistribute it and/or modify          #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 3 of the License, or       #
#    any later version.                                                      #
#                                                                            #
#    FemLab is distributed in the hope that it will be useful,               #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with FemLab.  If not, see <http://www.gnu.org/licenses/>.         #
##############################################################################

# This file includes the code for generating embedded line cells

export generate_embedded_cells!

# Remove joint cells from mesh and set up line cells as embedded cells
function generate_embedded_cells!(mesh::Mesh)

    newcells = []
    id = 0 
    for cell in mesh.cells
        if cell.shape.class==JOINT1D_SHAPE
            # link solid cell to line cells
            solid, line = cell.linked_cells
            line.linked_cells = [ solid ]
        else
            # save non joint1D cells
            id += 1
            cell.id = id
            push!(newcells, cell)
        end
    end

    # update mesh
    mesh.cells = newcells
    return mesh
end
