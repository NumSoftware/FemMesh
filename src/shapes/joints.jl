
# Conventional (two-layered) joint shapes: JLIN2, JLIN3, JLIN4, JTRI3, JQUAD4, JTRI6, JQUAD8
# ==========================================================================================

for (shape, npoints) in [ (:LIN2,4), (:LIN3,6), (:LIN4, 8), (:TRI3,6), (:QUAD4,8), (:TRI6,12), (:QUAD8,16) ] 
    # construction is based on modification of basic shapes
    name = "J$shape"
    makefunc = Symbol("Make$shape")
    j2shape = Symbol(name)

    ex = quote
        const $j2shape    = $makefunc()
        $j2shape.name     = $name
        $j2shape.family   = JOINT_SHAPE
        $j2shape.npoints  = $npoints
        $j2shape.basic_shape = $shape
        $j2shape.vtk_type    = VTK_POLY_VERTEX
        $j2shape.facet_idxs  = []
        $j2shape.edge_idxs   = []
        $j2shape.facet_shape = $shape
        export $j2shape
    end

    eval(ex)
end

# Three-layered joint shapes: J3LIN2, J3LIN3, J3LIN4, J3TRI3, J3QUAD4, J3TRI6, J3QUAD8
# ====================================================================================

for (shape, npoints) in [ (:LIN2,6), (:LIN3,9), (:LIN4, 8), (:TRI3,9), (:QUAD4,12), (:TRI6,18), (:QUAD8,24) ] 
    # construction is based on modification of basic shapes
    name = "J3$shape"
    makefunc = Symbol("Make$shape")
    j3shape = Symbol(name)

    ex = quote
        const $j3shape    = $makefunc()
        $j3shape.name     = $name
        $j3shape.family   = JOINT_SHAPE
        $j3shape.npoints  = $npoints
        $j3shape.basic_shape = $shape
        $j3shape.vtk_type    = VTK_POLY_VERTEX
        $j3shape.facet_idxs  = []
        $j3shape.edge_idxs   = []
        $j3shape.facet_shape = $shape
        export $j3shape
    end

    eval(ex)
end


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


