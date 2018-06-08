##############################################################################
#    FemLab - Finite Element Library                                         #
#    Copyright (C) Raul Durand <raul.durand at gmail.com>                    #
#    All rights reserved.                                                    #
#                                                                            #
#    This software is distributed WITHOUT ANY WARRANTY; without even         #
#    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR     #
#    PURPOSE. See the GNU General Public License for more details.           #
##############################################################################


# JLIN2 shape
# ==========

# constructor
function MakeJLIN2()
    shape             = ShapeType()
    shape.name        = "JLIN2"
    shape.family      = JOINT_SHAPE
    shape.ndim        = 1
    shape.npoints     = 2
    shape.basic_shape = LIN2
    shape.vtk_type    = VTK_POLY_VERTEX
    shape.facet_idxs  = []
    shape.edge_idxs   = []
    shape.facet_shape = LIN2
    shape.nat_coords  = coords_LIN2
    shape.quadrature  = Dict( 0 => LIN_IP2, 1 => ALL_IP1, 2 => LIN_IP2, 3 => LIN_IP3, 4 => LIN_IP4)
    shape.func        = shape_func_LIN2
    shape.deriv       = shape_deriv_LIN2
    return shape
end


# Registration
const JLIN2 = MakeJLIN2()
export JLIN2



# JLIN3 shape
# ===========

# constructor
function MakeJLIN3()
    shape             = ShapeType()
    shape.name        = "JLIN3"
    shape.family      = JOINT_SHAPE
    shape.ndim        = 1
    shape.npoints     = 3
    shape.basic_shape = LIN3
    shape.vtk_type    = VTK_POLY_VERTEX
    shape.facet_idxs  = []
    shape.edge_idxs   = []
    shape.facet_shape = LIN3
    shape.nat_coords  = coords_LIN3
    shape.quadrature  = Dict( 0 => LIN_IP2, 1 => ALL_IP1, 2 => LIN_IP2, 3 => LIN_IP3, 4 => LIN_IP4 )
    shape.func        = shape_func_LIN3
    shape.deriv       = shape_deriv_LIN3
    return shape
end


# Registration
const JLIN3 = MakeJLIN3()
export JLIN3



# JLIN4 shape
# ===========

# constructor
function MakeJLIN4()
    shape             = ShapeType()
    shape.name        = "JLIN4"
    shape.family      = JOINT_SHAPE
    shape.ndim        = 1
    shape.npoints     = 4
    shape.basic_shape = LIN4
    shape.vtk_type    = VTK_POLY_VERTEX
    shape.facet_idxs  = []
    shape.edge_idxs   = []
    shape.facet_shape = LIN4
    shape.nat_coords  = coords_LIN4
    shape.quadrature  = Dict( 0 => LIN_IP3, 1 => ALL_IP1, 2 => LIN_IP2, 3 => LIN_IP3, 4 => LIN_IP4 )
    shape.func        = shape_func_LIN4
    shape.deriv       = shape_deriv_LIN4
    return shape
end


# Registration
const JLIN4 = MakeJLIN4()
export JLIN4



# JTRI3 shape
# ===========

# constructor
function MakeJTRI3()
    shape             = ShapeType()
    shape.name        = "JTRI3"
    shape.family      = JOINT_SHAPE
    shape.ndim        = 2
    shape.npoints     = 3
    shape.basic_shape = TRI3
    shape.vtk_type    = VTK_POLY_VERTEX
    shape.facet_idxs  = []
    shape.edge_idxs   = []
    shape.facet_shape = TRI3
    shape.nat_coords  = coords_TRI3
    shape.quadrature  = Dict( 0 => TRI_IP3, 1 => ALL_IP1, 3 => TRI_IP3, 6 => TRI_IP6 )
    shape.func        = shape_func_TRI3
    shape.deriv       = shape_deriv_TRI3
    return shape
end

# Registration
const JTRI3 = MakeJTRI3()
export JTRI3




# JTRI6 shape
# ===========

# constructor
function MakeJTRI6()
    shape             = ShapeType()
    shape.name        = "JTRI6"
    shape.family      = JOINT_SHAPE
    shape.ndim        = 2
    shape.npoints     = 6
    shape.basic_shape = TRI3
    shape.vtk_type    = VTK_POLY_VERTEX
    shape.facet_idxs  = []
    shape.edge_idxs   = []
    shape.facet_shape = TRI6
    shape.nat_coords  = coords_TRI6
    shape.quadrature  = Dict( 0 => TRI_IP3, 1 => ALL_IP1, 3 => TRI_IP3, 6 => TRI_IP6 )
    shape.func        = shape_func_TRI6
    shape.deriv       = shape_deriv_TRI6
    return shape
end

# Registration
const JTRI6 = MakeJTRI6()
export JTRI6




# JQUAD4 shape
# ===========

# constructor
function MakeJQUAD4()
    shape             = ShapeType()
    shape.name        = "JQUAD4"
    shape.family      = JOINT_SHAPE
    shape.ndim        = 2
    shape.npoints     = 4
    shape.basic_shape = QUAD4
    shape.vtk_type    = VTK_POLY_VERTEX
    shape.facet_idxs  = []
    shape.edge_idxs   = []
    shape.facet_shape = QUAD4
    shape.nat_coords  = coords_QUAD4
    shape.quadrature  = Dict( 0 => QUAD_IP2, 1 => ALL_IP1, 4 => QUAD_IP2, 9 => QUAD_IP3 )
    shape.func        = shape_func_QUAD4
    shape.deriv       = shape_deriv_QUAD4
    return shape
end

# Registration
const JQUAD4 = MakeJQUAD4()
export JQUAD4




# JQUAD8 shape
# ===========

# constructor
function MakeJQUAD8()
    shape             = ShapeType()
    shape.name        = "JQUAD8"
    shape.family      = JOINT_SHAPE
    shape.ndim        = 2
    shape.npoints     = 8
    shape.basic_shape = QUAD8
    shape.vtk_type    = VTK_POLY_VERTEX
    shape.facet_idxs  = []
    shape.edge_idxs   = []
    shape.facet_shape = QUAD8
    shape.nat_coords  = coords_QUAD8
    shape.quadrature  = Dict( 0 => QUAD_IP3, 4 => QUAD_IP2, 9 => QUAD_IP3 )
    shape.func        = shape_func_QUAD8
    shape.deriv       = shape_deriv_QUAD8
    return shape
end

# Registration
const JQUAD8 = MakeJQUAD8()
export JQUAD8



# JLINK2 shape
# ===========

# constructor
function MakeJLINK2()
    shape             = ShapeType()
    shape.name        = "JLINK2"
    shape.family      = JOINT1D_SHAPE
    shape.ndim        = 1
    shape.npoints     = 2
    shape.basic_shape = LIN2
    shape.vtk_type    = VTK_POLY_VERTEX
    shape.facet_idxs  = []
    shape.edge_idxs   = []
    shape.facet_shape = ()
    shape.nat_coords  = coords_LIN2
    shape.quadrature  = Dict( 0 => LIN_IP2,  2 => LIN_IP2,  3 => LIN_IP3,   4 => LIN_IP4)
    shape.func        = shape_func_LIN2
    shape.deriv       = shape_deriv_LIN2
    return shape
end


# Registration
const JLINK2 = MakeJLINK2()
export JLINK2


# JLINK3 shape
# ===========

# constructor
function MakeJLINK3()
    shape             = ShapeType()
    shape.name        = "JLINK3"
    shape.family      = JOINT1D_SHAPE
    shape.ndim        = 1
    shape.npoints     = 3
    shape.basic_shape = LIN3
    shape.vtk_type    = VTK_POLY_VERTEX
    shape.facet_idxs  = []
    shape.edge_idxs   = []
    shape.facet_shape = ()
    shape.nat_coords  = coords_LIN3
    shape.quadrature  = Dict( 0 => LIN_IP3,  2 => LIN_IP2,  3 => LIN_IP3,   4 => LIN_IP4 )
    shape.func        = shape_func_LIN3
    shape.deriv       = shape_deriv_LIN3
    return shape
end


# Registration
const JLINK3 = MakeJLINK3()
export JLINK3


