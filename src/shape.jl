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


include("quadrature.jl")

export get_shape_tag, get_vtk_type, get_shape_from_vtk, get_ndim, get_name, get_basic_shape
export LIN2, LIN3, TRI3, TRI6, QUAD4, QUAD8, TET4, TET10, HEX8, HEX20, WED6, WED15
export LINK1, LINK2, LINK3, LIN4, TRI9, TRI10, QUAD9, QUAD12, QUAD16
export ShapeType
export shape_func, deriv_func, local_coords
export get_ip_coords
export is_solid, is_line, is_joint, is_joint1D, is_inside
export inverse_map, extrapolator
export coords_tri6, coords_tri9, coords_tri10, coords_quad4, coords_quad8, coords_wed6, coords_wed15 #...

@enum(ShapeType,
# Shapes compatible with VTK
POLYV = 2,
LIN2  = 3,
LIN3  = 21,
TRI3  = 5,
TRI6  = 22,
QUAD4 = 9,
QUAD8 = 23,
TET4  = 10,
TET10 = 24,
HEX8  = 12,
HEX20 = 25,
WED6  = 13,
# Shapes represented as polyvertex
LIN4   = 31,
TRI9   = 33,
TRI10  = 34,
QUAD9  = 37,
QUAD12 = 38,
QUAD16 = 39,
WED15  = 41,
# Embedded links
LINK1  = 51,
LINK2  = 52,
LINK3  = 53,
# Conventional joints
JLIN2  = 100 + 3,
JLIN3  = 100 + 21,
JLIN4  = 100 + 31,
JTRI3  = 100 + 5,
JTRI6  = 100 + 22,
JQUAD4 = 100 + 9,
JQUAD8 = 100 + 23,
)

function is_line(shape::ShapeType)
    shape in (LIN2, LIN3, LIN4) 
end

function is_joint1D(shape::ShapeType)
    Int(shape)>50 && Int(shape)<60
end

function is_joint(shape::ShapeType)
    Int(shape)>100
end

function is_solid(shape::ShapeType)
    Int(shape)<50 && !(shape in (LIN2, LIN3, LIN4, POLYV) )
end

function get_ndim(shape::ShapeType)
    # Returns the local dimension based on the shape geometry.
    # It does not match necessarily the space where the shape is used.
    
    if shape == LIN2  ; return 1 end
    if shape == LIN3  ; return 1 end
    if shape == LIN4  ; return 1 end
    if shape == LINK2 ; return 1 end
    if shape == LINK3 ; return 1 end
    if shape == TRI3  ; return 2 end
    if shape == TRI6  ; return 2 end
    if shape == TRI9  ; return 2 end
    if shape == TRI10 ; return 2 end
    if shape == QUAD4 ; return 2 end
    if shape == QUAD8 ; return 2 end
    if shape == QUAD9 ; return 2 end
    if shape == QUAD12; return 2 end
    if shape == QUAD16; return 2 end
    return 3
    error("Unknown shape ($shape)")
end

function get_nnodes(shape::ShapeType)
    # Returns the local dimension based on the shape geometry.
    # It does not match necessarily the space where the shape is used.
    
    if shape == LIN2  ; return  2 end
    if shape == LIN3  ; return  3 end
    if shape == LIN4  ; return  4 end
    if shape == LINK2 ; return  2 end
    if shape == LINK3 ; return  3 end
    if shape == TRI3  ; return  3 end
    if shape == TRI6  ; return  6 end
    if shape == TRI9  ; return  9 end
    if shape == TRI10 ; return 10 end
    if shape == QUAD4 ; return  4 end
    if shape == QUAD8 ; return  8 end
    if shape == QUAD9 ; return  9 end
    if shape == QUAD12; return 12 end
    if shape == QUAD16; return 16 end
    if shape == TET4  ; return  4 end
    if shape == TET10 ; return 10 end
    if shape == HEX8  ; return  8 end
    if shape == HEX20 ; return 20 end
    if shape == WED6  ; return  6 end
    if shape == WED15 ; return 15 end
    if shape>100;       return get_nnodes(Int(shape)-100)*2 end
    error("Unknown shape ($shape)")
end


#=
# Shapes compatible with VTK
const POLYV = 2
const LIN2  = 3
const LIN3  = 21
const TRI3  = 5
const TRI6  = 22
const QUAD4 = 9
const QUAD8 = 23
const TET4  = 10
const TET10 = 24
const HEX8  = 12
const HEX20 = 25
const WED6  = 13

# Shapes represented as polyvertex
const LIN4   = 31
const TRI9   = 33
const TRI10  = 34
const QUAD9  = 37
const QUAD12 = 38
const QUAD16 = 39
const WED15  = 41

# Embedded links
const LINK1  = 51
const LINK2  = 52
const LINK3  = 53

# Conventional joints
const JLIN2  = 100 + LIN2 
const JLIN3  = 100 + LIN3 
const JLIN4  = 100 + LIN4 
const JTRI3  = 100 + TRI3 
const JTRI6  = 100 + TRI6 
const JQUAD4 = 100 + QUAD4
const JQUAD8 = 100 + QUAD8

typealias ShapeType Int64
=#

SHAPE_TAG = Dict( LIN2  => 102, LIN3   => 103, LIN4 => 104, 
              TRI3  => 703, TRI6   => 706, TRI9 => 709, TRI10  => 710, 
              QUAD4 => 904, QUAD8  => 908, 
              TET4  => 704, TET10  => 710, HEX8  => 408, HEX20 => 420, 
              QUAD9 => 909, QUAD12 => 912, QUAD16 => 916, 
              LINK1 => 201, LINK2  => 202, LINK3  => 203 )

SHAPE_NAME = Dict( LIN2  => "LIN2" , LIN3   => "LIN3"  , LIN4   => "LIN4", 
               TRI3  => "TRI3" , TRI6   => "TRI6"  , TRI9   => "TRI9", TRI10 => "TRI10", 
               QUAD4 => "QUAD4", QUAD8  => "QUAD8" , QUAD9  => "QUAD9", QUAD12  => "QUAD12", QUAD16 => "QUAD16",
               TET4  => "TET4" , TET10  => "TET10" , HEX8   => "HEX8", HEX20 => "HEX20", 
               WED6  => "WED6" , WED15  => "WED15" , 
               QUAD9 => "QUAD9", QUAD12 => "QUAD12", QUAD16 => "QUAD16", 
               LINK1 => "LINK1", LINK2  => "LINK2" , LINK3  => "LINK3" )

function get_shape_tag(shape::ShapeType)
    return SHAPE_TAG[shape]
end

function get_name(shape::ShapeType)
    SHAPE_NAME[shape]
end

function get_basic_shape(shape::ShapeType)
    # Returs the basic shape of another
    
    if shape in [ LIN2, LIN3, LIN4 ] return LIN2 end
    if shape in [ TRI3, TRI6, TRI9, TRI10 ] return TRI3 end
    if shape in [ QUAD4, QUAD8, QUAD12, QUAD16 ] return QUAD4 end
    if shape in [ TET4, TET10 ] return TET4 end
    if shape in [ HEX8, HEX20 ] return HEX8 end
    if shape in [ WED6, WED15] return WED8 end

    if Int(shape) > 100 return ShapeType(Int(shape)-100) end # for joint shapes

    return shape # any other shape
end

function get_vtk_type(shape::ShapeType)::Int
    vtk_type = Int(shape)
    vtk_type < 50 && return vtk_type # Conventional
    return Int(POLYV)                # vtk_poly_vertex
end

function get_shape_from_vtk(vtk_type::Int64, npoints::Int64, ndim::Int64)::ShapeType

    if vtk_type!=2 return ShapeType(vtk_type) end # diff from poly_vertex

    if     npoints==1   return LINK1
    elseif npoints==2   return LINK2
    elseif npoints==3   return LINK3
    elseif npoints==4   return JLIN2
    elseif npoints==6   
        if ndim==2 return JLIN3 end
        if ndim==3 return JTRI3 end
    elseif npoints==8
        if ndim==2 return JLIN4  end
        if ndim==3 return JQUAD4 end
    elseif npoints==9   return QUAD9  # what about TRI9?
    elseif npoints==10  return TRI10
    elseif npoints==12  
        if ndim==2 return QUAD12 end
        if ndim==3 return JTRI6  end
    elseif npoints==15  
        if ndim==3 return WED15  end
    elseif npoints==16
        if ndim==2 return QUAD16 end
        if ndim==3 return JQUAD8 end
    elseif npoints==18  return JTRI9
    elseif npoints==20  return JTRI10
    end

    error("get_shape_from_vtk: Unknown shape for vtk_type $vtk_type and npoints $npoints with ndim $ndim")
    #return POLYV

end


coords_lin2 =
[ -1.0  1.0 
   1.0  1.0 ]

coords_lin3 =
[ -1.0  1.0 
   0.0  1.0 
   1.0  1.0 ]

coords_tri6 =
[ 0.0  0.0  1.0
  1.0  0.0  1.0
  0.0  1.0  1.0 
  0.5  0.0  1.0 
  0.5  0.5  1.0 
  0.0  0.5  1.0 ]


_1_3 = 1.0/3.0; _2_3 = 2.0/3.0

coords_tri9 =
[ 0.0   0.0  1.0
  1.0   0.0  1.0
  0.0   1.0  1.0
  
  _1_3   0.0  1.0
  _2_3  _1_3  1.0
   0.0  _2_3  1.0
  
  _2_3   0.0  1.0
  _1_3  _2_3  1.0
   0.0  _1_3  1.0 ]

coords_tri10 =
[ 0.0   0.0  1.0 
  1.0   0.0  1.0 
  0.0   1.0  1.0 
  
  _1_3   0.0  1.0 
  _2_3  _1_3  1.0 
   0.0  _2_3  1.0 
  
  _2_3   0.0  1.0 
  _1_3  _2_3  1.0 
   0.0  _1_3  1.0 
  
  _1_3  _1_3  1.0 ]


coords_quad4 = 
[ -1.0 -1.0  1.0
   1.0 -1.0  1.0
   1.0  1.0  1.0
  -1.0  1.0  1.0 ]

coords_quad8 = 
[ -1.0 -1.0  1.0
   1.0 -1.0  1.0
   1.0  1.0  1.0
  -1.0  1.0  1.0
   0.0 -1.0  1.0
   1.0  0.0  1.0
   0.0  1.0  1.0
  -1.0  0.0  1.0 ]

coords_hex8 = 
[ -1.0 -1.0 -1.0  1.0
   1.0 -1.0 -1.0  1.0
   1.0  1.0 -1.0  1.0
  -1.0  1.0 -1.0  1.0
  -1.0 -1.0  1.0  1.0
   1.0 -1.0  1.0  1.0
   1.0  1.0  1.0  1.0
  -1.0  1.0  1.0  1.0 ]

coords_hex20 = 
[ -1.0 -1.0 -1.0  1.0
   1.0 -1.0 -1.0  1.0
   1.0  1.0 -1.0  1.0
  -1.0  1.0 -1.0  1.0
  -1.0 -1.0  1.0  1.0
   1.0 -1.0  1.0  1.0
   1.0  1.0  1.0  1.0
  -1.0  1.0  1.0  1.0 
 
   0.0 -1.0 -1.0  1.0
   1.0  0.0 -1.0  1.0
   0.0  1.0 -1.0  1.0
  -0.0  0.0 -1.0  1.0
   0.0 -1.0  1.0  1.0
   1.0  0.0  1.0  1.0
   0.0  1.0  1.0  1.0
  -0.0  0.0  1.0  1.0

  -1.0 -1.0  0.0  1.0
   1.0 -1.0  0.0  1.0
   1.0  1.0  0.0  1.0
  -1.0  1.0  0.0  1.0
  ]

coords_tet4 = 
[  0.0  0.0  0.0  1.0 
   1.0  0.0  0.0  1.0 
   0.0  1.0  0.0  1.0 
   0.0  0.0  1.0  1.0 ]

coords_tet10 = 
[  0.0  0.0  0.0  1.0 
   1.0  0.0  0.0  1.0 
   0.0  1.0  0.0  1.0 
   0.0  0.0  1.0  1.0 

   0.5  0.0  0.0  1.0 
   0.5  0.5  0.0  1.0 
   0.0  0.5  0.0  1.0 
   0.0  0.0  0.5  1.0 
   0.5  0.0  0.5  1.0 
   0.0  0.5  0.5  1.0 ]

coords_wed6 = 
[ 0.0  0.0 -1.0  1.0
  1.0  0.0 -1.0  1.0
  0.0  1.0 -1.0  1.0
  0.0  0.0  1.0  1.0
  1.0  0.0  1.0  1.0
  0.0  1.0  1.0  1.0 ]

coords_wed15 = 
[ 0.0  0.0 -1.0  1.0
  1.0  0.0 -1.0  1.0
  0.0  1.0 -1.0  1.0
  0.0  0.0  1.0  1.0
  1.0  0.0  1.0  1.0
  0.0  1.0  1.0  1.0

  0.5  0.0 -1.0  1.0
  0.5  0.5 -1.0  1.0
  0.0  0.5 -1.0  1.0
  0.5  0.0  1.0  1.0
  0.5  0.5  1.0  1.0
  0.0  0.5  1.0  1.0
  0.0  0.0  0.0  1.0
  1.0  0.0  0.0  1.0
  0.0  1.0  0.0  1.0 ]

function get_local_coords(st::ShapeType)
    if     st == LIN2   return coords_lin2
    elseif st == LIN3   return coords_lin3
    elseif st == TRI6   return coords_tri6
    elseif st == QUAD4  return coords_quad4
    elseif st == TET10  return coords_tet10
    elseif st == HEX20  return coords_hex20
    elseif st == WED6   return coords_wed6
    elseif st == WED15  return coords_wed15
    end

    error("get_local_coords: No local coordinates implemented for shape $st")
    return 0
end


#immutable Typed{N} end

shape_func(shape::ShapeType, R::Array{Float64,1}) = shape_func(Val{shape}, R)
deriv_func(shape::ShapeType, R::Array{Float64,1}) = deriv_func(Val{shape}, R)

function shape_func(::Type{Val{LIN2}}, R::Array{Float64,1})
    r = R[1]
    N = Array(Float64, 2)
    N[1] = 0.5*(1-r)
    N[2] = 0.5*(1+r)
    return N
end

function deriv_func(::Type{Val{LIN2}}, R::Array{Float64,1})
    D = Array(Float64, 1, 2)
    D[1,1] = -0.5
    D[1,2] =  0.5
    return D
end

function shape_func(::Type{Val{LIN3}}, R::Array{Float64,1})
    r = R[1]
    N = Array(Float64, 3)
    N[1] = 0.5*(r*r - r)
    N[2] = 0.5*(r*r + r)
    N[3] = 1.0 - r*r
    return N
end

function deriv_func(::Type{Val{LIN3}}, R::Array{Float64,1})
    r = R[1]
    D = Array(Float64, 1, 3)
    D[1, 1] = r - 0.5
    D[1, 2] = r + 0.5
    D[1, 3] = -2.0*r
    return D
end

function shape_func(::Type{Val{LIN4}}, R::Array{Float64,1})
    #   (-1)            '   (+1)
    #    @------@-----@------@  --> r
    #    1      3     4      2

    r = R[1]
    N = Array(Float64, 4)
    N[1] = 1./16.*( -9.*r^3 + 9.*r*r +     r - 1.)
    N[2] = 1./16.*(  9.*r^3 + 9.*r*r -     r - 1.)
    N[3] = 1./16.*( 27.*r^3 - 9.*r*r - 27.*r + 9.)
    N[4] = 1./16.*(-27.*r^3 - 9.*r*r + 27.*r + 9.)
    return N
end

function deriv_func(::Type{Val{LIN4}}, R::Array{Float64,1})
    r = R[1]
    D = Array(Float64, 1, 4)
    D[1,1] = 1./16.*( -27.*r*r + 18.*r + 1. )
    D[1,2] = 1./16.*(  27.*r*r + 18.*r - 1. )
    D[1,3] = 1./16.*(  81.*r*r - 18.*r - 27.)
    D[1,4] = 1./16.*( -81.*r*r - 18.*r + 27.)
    return D
end

function shape_func(::Type{Val{TRI3}}, R::Array{Float64,1})
    #    s
    #    ^
    #    |
    #  3
    #    @,(0,1)
    #    | ',
    #    |   ',
    #    |     ',
    #    |       ',
    #    |         ',
    #    |           ',
    #    |             ',
    #    |               ',
    #    |(0,0)            ', (1,0)
    #    @-------------------@  --> r
    #  1                      2
    #
    r, s = R[1:2]
    N = Array(Float64,3)
    N[1] = 1.0-r-s
    N[2] = r
    N[3] = s
    return N
end

function deriv_func(::Type{Val{TRI3}}, R::Array{Float64,1})
    r, s = R
    D = Array(Float64, 2, 3)
    D[1,1] = -1.0;    D[2,1] = -1.0
    D[1,2] =  1.0;    D[2,2] =  0.0
    D[1,3] =  0.0;    D[2,3] =  1.0
    return D
end

function shape_func(::Type{Val{TRI6}}, R::Array{Float64,1})
    #    s


    #    ^
    #    |
    #  2
    #    @,(0,1)
    #    | ',
    #    |   ',
    #    |     ',
    #    |       ',  4
    #  5 @         '@
    #    |           ',
    #    |             ',
    #    |               ',
    #    |(0,0)            ', (1,0)
    #    @---------@---------@  --> r
    #  0           3          1
    #
    r, s = R

    N = Array(Float64,6)
    N[1] = 1.0-(r+s)*(3.0-2.0*(r+s))
    N[2] = r*(2.0*r-1.0)
    N[3] = s*(2.0*s-1.0)
    N[4] = 4.0*r*(1.0-(r+s))
    N[5] = 4.0*r*s
    N[6] = 4.0*s*(1.0-(r+s))

    return N
end

function deriv_func(::Type{Val{TRI6}}, R::Array{Float64,1})
    r, s = R

    D = Array(Float64, 2, 6)
    D[1,1] = -3.0 + 4.0 * (r + s);       D[2,1] = -3.0 + 4.0*(r + s)
    D[1,2] =  4.0 * r - 1.;              D[2,2] =  0.0
    D[1,3] =  0.0;                       D[2,3] =  4.0 * s - 1.0
    D[1,4] =  4.0 - 8.0 * r - 4.0 * s;   D[2,4] = -4.0 * r
    D[1,5] =  4.0 * s;                   D[2,5] =  4.0 * r
    D[1,6] = -4.0 * s;                   D[2,6] =  4.0 - 4.0 * r - 8.0*s

    return D
end

function shape_func(::Type{Val{QUAD4}}, R::Array{Float64,1})
    #     4                        3
    #       @--------------------@
    #       |               (1,1)|
    #       |       s ^          |
    #       |         |          |
    #       |         |          |
    #       |         +----> r   |
    #       |       (0,0)        |
    #       |                    |
    #       |                    |
    #       |(-1,-1)             |
    #       @--------------------@
    #     1                        2
    #
    r, s = R[1:2]
    N = Array(Float64,4)
    N[1] = 0.25*(1.0-r-s+r*s)
    N[2] = 0.25*(1.0+r-s-r*s)
    N[3] = 0.25*(1.0+r+s+r*s)
    N[4] = 0.25*(1.0-r+s-r*s)
    return N
end

function deriv_func(::Type{Val{QUAD4}}, R::Array{Float64,1})
    r, s = R[1:2]
    D = Array(Float64, 2, 4)
    D[1,1] = 0.25*(-1.0+s);   D[2,1] = 0.25*(-1.0+r)
    D[1,2] = 0.25*(+1.0-s);   D[2,2] = 0.25*(-1.0-r)
    D[1,3] = 0.25*(+1.0+s);   D[2,3] = 0.25*(+1.0+r)
    D[1,4] = 0.25*(-1.0-s);   D[2,4] = 0.25*(+1.0-r)
    return D
end

function shape_func(::Type{Val{QUAD8}}, R::Array{Float64,1})
    #     4           7            3
    #       @---------@----------@
    #       |               (1,1)|
    #       |       s ^          |
    #       |         |          |
    #       |         |          |
    #     8 @         +----> r   @ 6
    #       |       (0,0)        |
    #       |                    |
    #       |                    |
    #       |(-1,-1)             |
    #       @---------@----------@
    #     1           5            2
    #
    r, s = R[1:2]
    N = Array(Float64,8)
    rp1=1.0+r; rm1=1.0-r;
    sp1=1.0+s; sm1=1.0-s;
    N[1] = 0.25*rm1*sm1*(rm1+sm1-3.0)
    N[2] = 0.25*rp1*sm1*(rp1+sm1-3.0)
    N[3] = 0.25*rp1*sp1*(rp1+sp1-3.0)
    N[4] = 0.25*rm1*sp1*(rm1+sp1-3.0)
    N[5] = 0.50*sm1*(1.0-r*r)
    N[6] = 0.50*rp1*(1.0-s*s)
    N[7] = 0.50*sp1*(1.0-r*r)
    N[8] = 0.50*rm1*(1.0-s*s)
    return N
end

function deriv_func(::Type{Val{QUAD8}}, R::Array{Float64,1})
    r, s = R[1:2]
    D = Array(Float64, 2, 8)
    rp1=1.0+r; rm1=1.0-r
    sp1=1.0+s; sm1=1.0-s

    D[1,1] = -0.25 * sm1 * (rm1 + rm1 + sm1 - 3.0)
    D[1,2] =  0.25 * sm1 * (rp1 + rp1 + sm1 - 3.0)
    D[1,3] =  0.25 * sp1 * (rp1 + rp1 + sp1 - 3.0)
    D[1,4] = -0.25 * sp1 * (rm1 + rm1 + sp1 - 3.0)
    D[1,5] = -r * sm1
    D[1,6] =  0.50 * (1.0 - s * s)
    D[1,7] = -r * sp1
    D[1,8] = -0.5 * (1.0 - s * s)

    D[2,1] = -0.25 * rm1 * (sm1 + rm1 + sm1 - 3.0)
    D[2,2] = -0.25 * rp1 * (sm1 + rp1 + sm1 - 3.0)
    D[2,3] =  0.25 * rp1 * (sp1 + rp1 + sp1 - 3.0)
    D[2,4] =  0.25 * rm1 * (sp1 + rm1 + sp1 - 3.0)
    D[2,5] = -0.50 * (1.0 - r * r)
    D[2,6] = -s * rp1
    D[2,7] =  0.50 * (1.0 - r * r)
    D[2,8] = -s * rm1
    return D
end

function shape_func(::Type{Val{QUAD9}}, R::Array{Float64,1})
    #     4           7            3
    #       @---------@----------@
    #       |               (1,1)|
    #       |       s ^          |
    #       |         |          |
    #       |         | 9        |
    #     8 @         +----> r   @ 6
    #       |       (0,0)        |
    #       |                    |
    #       |                    |
    #       |(-1,-1)             |
    #       @---------@----------@
    #     1           5            2
    #
    r, s = R[1:2]
    N = Array(Float64,9)
    rp1=1.0+r; rm1=1.0-r;
    sp1=1.0+s; sm1=1.0-s;
    N[1] =  0.25*r*s*sm1*rm1
    N[2] =  0.25*r*s*sm1*rp1
    N[3] =  0.25*r*s*sp1*rp1
    N[4] =  0.25*r*s*sp1*rm1
    N[5] = -0.50*s*sm1*rp1*rm1
    N[6] = -0.50*r*sp1*sm1*rp1
    N[7] = -0.50*s*sp1*rp1*rm1
    N[8] = -0.50*r*sp1*sm1*rm1
    N[9] = sp1*sm1*rp1*rm1
    return N
end

function deriv_func(::Type{Val{QUAD9}}, R::Array{Float64,1})
    r, s = R[1:2]
    D = Array(Float64, 2,9)
    rp1=1.0+r; rm1=1.0-r
    sp1=1.0+s; sm1=1.0-s

    #RP = 1. + r
    #RM = 1. - r
    #SP = 1. + s
    #SM = 1. - s
	D[1,1] = (r + r - 1.0) * s * (s - 1.0) / 4.0
	D[1,2] = (r + r + 1.0) * s * (s - 1.0) / 4.0
	D[1,3] = (r + r + 1.0) * s * (s + 1.0) / 4.0
	D[1,4] = (r + r - 1.0) * s * (s + 1.0) / 4.0
	D[1,5] = -(r + r) * s * (s - 1.0) / 2.0
	D[1,6] = -(r + r + 1.0) * (s*s - 1.0) / 2.0
	D[1,7] = -(r + r) * s * (s + 1.0) / 2.0
	D[1,8] = -(r + r - 1.0) * (s*s - 1.0) / 2.0
	D[1,9] = 2.0 * r * (s*s - 1.0)

	D[2,1] = r * (r - 1.0) * (s + s - 1.0) / 4.0
	D[2,2] = r * (r + 1.0) * (s + s - 1.0) / 4.0
	D[2,3] = r * (r + 1.0) * (s + s + 1.0) / 4.0
	D[2,4] = r * (r - 1.0) * (s + s + 1.0) / 4.0
	D[2,5] = -(r*r - 1.0) * (s + s - 1.0) / 2.0
	D[2,6] = -r * (r + 1.0) * (s + s) / 2.0
	D[2,7] = -(r*r - 1.0) * (s + s + 1.0) / 2.0
	D[2,8] = -r * (r - 1.0) * (s + s) / 2.0
	D[2,9] = 2.0 * s * (r*r - 1.0)

	#D[1,1] = (r + rm1) * s * sm1 / 4.0
	#D[1,2] = (r + rp1) * s * sm1 / 4.0
	#D[1,3] = (r + rp1) * s * sp1 / 4.0
	#D[1,4] = (r + rm1) * s * sp1 / 4.0
	#D[1,5] = -(r + r) * s * (sm1) / 2.0
	#D[1,6] = -(r + rp1) * (s*s - 1.0) / 2.0
	#D[1,7] = -(r + r) * s * (sp1) / 2.0
	#D[1,8] = -(r + rm1) * (s*s - 1.0) / 2.0
	#D[1,9] = 2.0 * r * (s*s - 1.0)
#
	#D[2,1] = r * rm1 * (s + sm1) / 4.0
	#D[2,2] = r * rp1 * (s + sm1) / 4.0
	#D[2,3] = r * rp1 * (s + sp1) / 4.0
	#D[2,4] = r * rm1 * (s + sp1) / 4.0
	#D[2,5] = -(r*r - 1.0) * (s + sm1) / 2.0
	#D[2,6] = -r * (rp1) * (s + s) / 2.0
	#D[2,7] = -(r*r - 1.0) * (s + sp1) / 2.0
	#D[2,8] = -r * (rm1) * (s + s) / 2.0
	#D[2,9] = 2.0 * s * (r*r - 1.0)

    return D
end


function shape_func(::Type{Val{QUAD12}}, R::Array{Float64,1})
    #  
    #    4      11       7        3
    #      @-----@-------@------@
    #      |               (1,1)|
    #      |       s ^          |
    #    8 @         |          @ 10
    #      |         |          |
    #      |         +----> r   |
    #      |       (0,0)        |
    #   12 @                    @ 6
    #      |                    |
    #      |(-1,-1)             |
    #      @-----@-------@------@
    #    1       5       9        2
    #  
    r, s = R[1:2]
    N = Array(Float64,12)

    RM = 1. - r
    RP = 1. + r
    SM = 1. - s
    SP = 1. + s
    N[1]  = RM*SM*( 9.*(r*r + s*s) - 10.)/32.
    N[2]  = RP*SM*( 9.*(r*r + s*s) - 10.)/32.
    N[3]  = RP*SP*( 9.*(r*r + s*s) - 10.)/32.
    N[4]  = RM*SP*( 9.*(r*r + s*s) - 10.)/32.
    N[5]  = 9.*(1. - r*r)*(1. - 3.*r)*SM/32.
    N[6]  = 9.*(1. - s*s)*(1. - 3.*s)*RP/32.
    N[7]  = 9.*(1. - r*r)*(1. + 3.*r)*SP/32.
    N[8]  = 9.*(1. - s*s)*(1. + 3.*s)*RM/32.
    N[9]  = 9.*(1. - r*r)*(1. + 3.*r)*SM/32.
    N[10] = 9.*(1. - s*s)*(1. + 3.*s)*RP/32.
    N[11] = 9.*(1. - r*r)*(1. - 3.*r)*SP/32.
    N[12] = 9.*(1. - s*s)*(1. - 3.*s)*RM/32.

    return N
end

function deriv_func(::Type{Val{QUAD12}}, R::Array{Float64,1})
    r, s = R[1:2]
    D = Array(Float64, 2, 12)

    RP = 1. + r
    RM = 1. - r
    SP = 1. + s
    SM = 1. - s

    D[1,1]  =  SM*(9.*(2.*r - 3.*r*r - s*s) + 10.)/32.
    D[1,2]  =  SM*(9.*(2.*r + 3.*r*r + s*s) - 10.)/32.
    D[1,3]  =  SP*(9.*(2.*r + 3.*r*r + s*s) - 10.)/32.
    D[1,4]  =  SP*(9.*(2.*r - 3.*r*r - s*s) + 10.)/32.
    D[1,5]  =  9.*SM*(9.*r*r - 2.*r - 3.)/32.
    D[1,6]  =  9.*(1. - s*s)*(1. - 3.*s)/32.
    D[1,7]  =  9.*SP*(-9.*r*r - 2.*r + 3.)/32.
    D[1,8]  = -9.*(1. - s*s)*(1. + 3.*s)/32.
    D[1,9]  =  9.*SM*(-9.*r*r - 2.*r + 3.)/32.
    D[1,10] =  9.*(1. - s*s)*(1. + 3.*s)/32.
    D[1,11] =  9.*SP*(9.*r*r - 2.*r - 3.)/32.
    D[1,12] = -9.*(1. - s*s)*(1. - 3.*s)/32.
    D[2,1]  =  RM*(9.*(2.*s - 3.*s*s - r*r) + 10.)/32.
    D[2,2]  =  RP*(9.*(2.*s - 3.*s*s - r*r) + 10.)/32.
    D[2,3]  =  RP*(9.*(2.*s + 3.*s*s + r*r) - 10.)/32.
    D[2,4]  =  RM*(9.*(2.*s + 3.*s*s + r*r) - 10.)/32.
    D[2,5]  = -9.*(1. - r*r)*(1. - 3.*r)/32.
    D[2,6]  =  9.*RP*(9.*s*s - 2.*s - 3.)/32.
    D[2,7]  =  9.*(1. - r*r)*(1. + 3.*r)/32.
    D[2,8]  =  9.*RM*(-9.*s*s - 2.*s + 3.)/32.
    D[2,9]  = -9.*(1. - r*r)*(1. + 3.*r)/32.
    D[2,10] =  9.*RP*(-9.*s*s - 2.*s + 3.)/32.
    D[2,11] =  9.*(1. - r*r)*(1. - 3.*r)/32.
    D[2,12] =  9.*RM*(9.*s*s - 2.*s - 3.)/32.

    return D
end

function shape_func(::Type{Val{HEX8}}, R::Array{Float64,1})
    # Local IDs
    #                  Nodes                                   Faces
    #     z
    #     |           5                  8
    #    ,+--y         @________________@                    +________________+
    #  x'            ,'|              ,'|                  ,'|              ,'|
    #              ,'  |            ,'  |                ,'  |  ___       ,'  |
    #            ,'    |          ,'    |              ,'    |,'5,'  [0],'    |
    #      6   ,'      |      7 ,'      |            ,'      |~~~     ,'      |
    #        @'===============@'        |          +'===============+'  ,'|   |
    #        |         |      |         |          |   ,'|   |      |   |3|   |
    #        |         |      |         |          |   |2|   |      |   |,'   |
    #        |       1 @______|_________@          |   |,'   +______|_________+
    #        |       ,'       |       ,' 4         |       ,'       |       ,'
    #        |     ,'         |     ,'             |     ,' [1]  ___|     ,'
    #        |   ,'           |   ,'               |   ,'      ,'4,'|   ,'
    #        | ,'             | ,'                 | ,'        ~~~  | ,'
    #        @________________@'                   +________________+'
    #      2                   3

    r, s, t = R[1:3]
    N = Array(Float64,8)
    N[1] = 0.125*(1.0-r-s+r*s-t+s*t+r*t-r*s*t)
    N[2] = 0.125*(1.0+r-s-r*s-t+s*t-r*t+r*s*t)
    N[3] = 0.125*(1.0+r+s+r*s-t-s*t-r*t-r*s*t)
    N[4] = 0.125*(1.0-r+s-r*s-t-s*t+r*t+r*s*t)
    N[5] = 0.125*(1.0-r-s+r*s+t-s*t-r*t+r*s*t)
    N[6] = 0.125*(1.0+r-s-r*s+t-s*t+r*t-r*s*t)
    N[7] = 0.125*(1.0+r+s+r*s+t+s*t+r*t+r*s*t)
    N[8] = 0.125*(1.0-r+s-r*s+t+s*t-r*t-r*s*t)
    return N
end

function deriv_func(::Type{Val{HEX8}}, R::Array{Float64,1})
    r, s, t = R
    st = s*t
    rt = r*t
    rs = r*s
    D = Array(Float64, 3, 8)
    D[1,1] = -1.0+s+t-st;   D[2,1]=-1.0+r+t-rt;   D[3,1]=-1.0+r+s-rs
    D[1,2] = +1.0-s-t+st;   D[2,2]=-1.0-r+t+rt;   D[3,2]=-1.0-r+s+rs
    D[1,3] = +1.0+s-t-st;   D[2,3]=+1.0+r-t-rt;   D[3,3]=-1.0-r-s-rs
    D[1,4] = -1.0-s+t+st;   D[2,4]=+1.0-r-t+rt;   D[3,4]=-1.0+r-s+rs
    D[1,5] = -1.0+s-t+st;   D[2,5]=-1.0+r-t+rt;   D[3,5]=+1.0-r-s+rs
    D[1,6] = +1.0-s+t-st;   D[2,6]=-1.0-r-t-rt;   D[3,6]=+1.0+r-s-rs
    D[1,7] = +1.0+s+t+st;   D[2,7]=+1.0+r+t+rt;   D[3,7]=+1.0+r+s+rs
    D[1,8] = -1.0-s-t-st;   D[2,8]=+1.0-r+t-rt;   D[3,8]=+1.0-r+s-rs
    D = 0.125*D

    return D
end

function shape_func(::Type{Val{HEX20}}, R::Array{Float64,1})
    # Local IDs
    #                   Vertices                               Faces
    #     t
    #     |           5        16        8
    #    ,+--s         @-------@--------@                   +----------------+
    #  r'            ,'|              ,'|                 ,'|              ,'|
    #           13 @'  |         15 ,'  |               ,'  |  ___       ,'  |
    #            ,'    |17        ,@    |20           ,'    |,'6,'  [1],'    |
    #      6   ,'      @      7 ,'      @           ,'      |~~~     ,'      |
    #        @'=======@=======@'        |         +'===============+'  ,'|   |
    #        |      14 |      |         |         |   ,'|   |      |   |4|   |
    #        |         |      |  12     |         |   |3|   |      |   |,'   |
    #     18 |       1 @- - - | @- - - -@         |   |,'   +- - - | +- - - -+
    #        @       ,'       @       ,' 4        |       ,'       |       ,'
    #        |   9 @'      19 |     ,'            |     ,' [2]  ___|     ,'
    #        |   ,'           |   ,@ 11           |   ,'      ,'5,'|   ,'
    #        | ,'             | ,'                | ,'        ~~~  | ,'
    #        @-------@--------@'                  +----------------+'
    #      2        10         3

    r, s, t = R

    rp1=1.0+r; rm1=1.0-r
    sp1=1.0+s; sm1=1.0-s
    tp1=1.0+t; tm1=1.0-t

    N = Array(Float64,20)
    N[ 1] = 0.125*rm1*sm1*tm1*(-r-s-t-2.0)
    N[ 2] = 0.125*rp1*sm1*tm1*( r-s-t-2.0)
    N[ 3] = 0.125*rp1*sp1*tm1*( r+s-t-2.0)
    N[ 4] = 0.125*rm1*sp1*tm1*(-r+s-t-2.0)
    N[ 5] = 0.125*rm1*sm1*tp1*(-r-s+t-2.0)
    N[ 6] = 0.125*rp1*sm1*tp1*( r-s+t-2.0)
    N[ 7] = 0.125*rp1*sp1*tp1*( r+s+t-2.0)
    N[ 8] = 0.125*rm1*sp1*tp1*(-r+s+t-2.0)
    N[ 9] = 0.25*(1.0-r*r)*sm1*tm1
    N[10] = 0.25*rp1*(1.0-s*s)*tm1
    N[11] = 0.25*(1.0-r*r)*sp1*tm1
    N[12] = 0.25*rm1*(1.0-s*s)*tm1
    N[13] = 0.25*(1.0-r*r)*sm1*tp1
    N[14] = 0.25*rp1*(1.0-s*s)*tp1
    N[15] = 0.25*(1.0-r*r)*sp1*tp1
    N[16] = 0.25*rm1*(1.0-s*s)*tp1
    N[17] = 0.25*rm1*sm1*(1.0-t*t)
    N[18] = 0.25*rp1*sm1*(1.0-t*t)
    N[19] = 0.25*rp1*sp1*(1.0-t*t)
    N[20] = 0.25*rm1*sp1*(1.0-t*t)
    return N
end

function deriv_func(::Type{Val{HEX20}}, R::Array{Float64,1})
    r, s, t = R

    rp1=1.0+r; rm1=1.0-r
    sp1=1.0+s; sm1=1.0-s
    tp1=1.0+t; tm1=1.0-t

    D = Array(Float64, 3, 20)
    # Derivatives with respect to r
    D[1, 1] = -0.125*sm1*tm1*(-r-s-t-2)-0.125*rm1*sm1*tm1
    D[1, 2] =  0.125*sm1*tm1*( r-s-t-2)+0.125*rp1*sm1*tm1
    D[1, 3] =  0.125*sp1*tm1*( r+s-t-2)+0.125*rp1*sp1*tm1
    D[1, 4] = -0.125*sp1*tm1*(-r+s-t-2)-0.125*rm1*sp1*tm1
    D[1, 5] = -0.125*sm1*tp1*(-r-s+t-2)-0.125*rm1*sm1*tp1
    D[1, 6] =  0.125*sm1*tp1*( r-s+t-2)+0.125*rp1*sm1*tp1
    D[1, 7] =  0.125*sp1*tp1*( r+s+t-2)+0.125*rp1*sp1*tp1
    D[1, 8] = -0.125*sp1*tp1*(-r+s+t-2)-0.125*rm1*sp1*tp1
    D[1, 9] = -0.5*r*sm1*tm1
    D[1,10] =  0.25*(1-s*s)*tm1
    D[1,11] = -0.5*r*sp1*tm1
    D[1,12] = -0.25*(1-s*s)*tm1
    D[1,13] = -0.5*r*sm1*tp1
    D[1,14] =  0.25*(1-s*s)*tp1
    D[1,15] = -0.5*r*sp1  *tp1
    D[1,16] = -0.25*(1-s*s)*tp1
    D[1,17] = -0.25*sm1*(1-t*t)
    D[1,18] =  0.25*sm1*(1-t*t)
    D[1,19] =  0.25*sp1*(1-t*t)
    D[1,20] = -0.25*sp1*(1-t*t)

    # Derivatives with respect to s
    D[2, 1] = -0.125*rm1*tm1*(-r-s-t-2)-0.125*rm1*sm1*tm1
    D[2, 2] = -0.125*rp1*tm1*( r-s-t-2)-0.125*rp1*sm1*tm1
    D[2, 3] =  0.125*rp1*tm1*( r+s-t-2)+0.125*rp1*sp1*tm1
    D[2, 4] =  0.125*rm1*tm1*(-r+s-t-2)+0.125*rm1*sp1*tm1
    D[2, 5] = -0.125*rm1*tp1*(-r-s+t-2)-0.125*rm1*sm1*tp1
    D[2, 6] = -0.125*rp1*tp1*( r-s+t-2)-0.125*rp1*sm1*tp1
    D[2, 7] =  0.125*rp1*tp1*( r+s+t-2)+0.125*rp1*sp1*tp1
    D[2, 8] =  0.125*rm1*tp1*(-r+s+t-2)+0.125*rm1*sp1*tp1
    D[2, 9] = -0.25*(1-r*r)*tm1
    D[2,10] = -0.5*s*rp1*tm1
    D[2,11] =  0.25*(1-r*r)*tm1
    D[2,12] = -0.5*s*rm1*tm1
    D[2,13] = -0.25*(1-r*r)*tp1
    D[2,14] = -0.5*s*rp1*tp1
    D[2,15] =  0.25*(1-r*r)*tp1
    D[2,16] = -0.5*s*rm1*tp1
    D[2,17] = -0.25*rm1*(1-t*t)
    D[2,18] = -0.25*rp1*(1-t*t)
    D[2,19] =  0.25*rp1*(1-t*t)
    D[2,20] =  0.25*rm1*(1-t*t)

    # Derivatives with respect to t
    D[3, 1] = -0.125*rm1*sm1*(-r-s-t-2)-0.125*rm1*sm1*tm1
    D[3, 2] = -0.125*rp1*sm1*( r-s-t-2)-0.125*rp1*sm1*tm1
    D[3, 3] = -0.125*rp1*sp1*( r+s-t-2)-0.125*rp1*sp1*tm1
    D[3, 4] = -0.125*rm1*sp1*(-r+s-t-2)-0.125*rm1*sp1*tm1
    D[3, 5] =  0.125*rm1*sm1*(-r-s+t-2)+0.125*rm1*sm1*tp1
    D[3, 6] =  0.125*rp1*sm1*( r-s+t-2)+0.125*rp1*sm1*tp1
    D[3, 7] =  0.125*rp1*sp1*( r+s+t-2)+0.125*rp1*sp1*tp1
    D[3, 8] =  0.125*rm1*sp1*(-r+s+t-2)+0.125*rm1*sp1*tp1
    D[3, 9] = -0.25*(1-r*r)*sm1
    D[3,10] = -0.25*rp1*(1-s*s)
    D[3,11] = -0.25*(1-r*r)*sp1
    D[3,12] = -0.25*rm1*(1-s*s)
    D[3,13] =  0.25*(1-r*r)*sm1
    D[3,14] =  0.25*rp1*(1-s*s)
    D[3,15] =  0.25*(1-r*r)*sp1
    D[3,16] =  0.25*rm1*(1-s*s)
    D[3,17] = -0.5*t*rm1*sm1
    D[3,18] = -0.5*t*rp1*sm1
    D[3,19] = -0.5*t*rp1*sp1
    D[3,20] = -0.5*t*rm1*sp1

    return D
end

function shape_func(::Type{Val{TET4}}, R::Array{Float64,1})
    r, s, t = R

    N = Array(Float64,4)
    N[1] = 1.0-r-s-t
    N[2] = r
    N[3] = s
    N[4] = t
    return N
end

function deriv_func(::Type{Val{TET4}}, R::Array{Float64,1})
    r, s, t = R

    D = Array(Float64, 3, 4)
    D[1,1] = -1.0;   D[2,1]= -1.0;   D[3,1]= -1.0
    D[1,2] =  1.0;   D[2,2]=  0.0;   D[3,2]=  0.0
    D[1,3] =  0.0;   D[2,3]=  1.0;   D[3,3]=  0.0
    D[1,4] =  0.0;   D[2,4]=  0.0;   D[3,4]=  1.0
    return D
end

function shape_func(::Type{Val{TET10}}, R::Array{Float64,1})

    #                       t
    #                       |
    #                       |
    #                       | 4
    #                       @,
    #                      /|`
    #                      ||  `,
    #                     / |    ',
    #                     | |      \
    #                    /  |       `.
    #                    |  |         `,  10
    #                   /   @ 8         `@
    #                   |   |             \
    #                  /    |              `.
    #                  |    |                ',
    #               9 @     |                  \
    #                 |     @.,,_       7       `.
    #                |     / 1   ``'-.,,@_        `.
    #                |    /              ``''-.,,_  ', 3
    #               |    /                        ``'@.,,,
    #               |   '                       ,.-``     ``''- s
    #              |  ,@ 5                 _,-'`
    #              ' /                 ,.'`
    #             | /             _.@``
    #             '/          ,-'`   6
    #            |/      ,.-``
    #            /  _,-``
    #          .@ '`
    #         / 2
    #        /
    #       /
    #      r
    # 

    r, s, t = R

    N = Array(Float64,10)

    u = 1.0 - r - s - t

    # corners
    N[1] = u*(2.0*u - 1.0)
    N[2] = r*(2.0*r - 1.0)
    N[3] = s*(2.0*s - 1.0)
    N[4] = t*(2.0*t - 1.0)

    # midedge
    N[5] = 4.0 * u * r
    N[6] = 4.0 * r * s
    N[7] = 4.0 * s * u
    N[8] = 4.0 * u * t
    N[9] = 4.0 * r * t
    N[10] = 4.0 * s * t

    return N
end

function deriv_func(::Type{Val{TET10}}, R::Array{Float64,1})
    r, s, t = R

    D = Array(Float64, 3, 10)

    # r-derivatives: dN0/dr to dN9/dr
    D[1,1]  =  4.0*(r + s + t) - 3.0
    D[1,2]  =  4.0*r - 1.0
    D[1,3]  =  0.0
    D[1,4]  =  0.0
    D[1,5]  =  4.0 - 8.0*r - 4.0*s - 4.0*t
    D[1,6]  =  4.0*s
    D[1,7]  = -4.0*s
    D[1,8]  = -4.0*t
    D[1,9]  =  4.0*t
    D[1,10] =  0.0

    # s-derivatives: dN0/ds to dN9/ds
    D[2,1]  =  4.0*(r + s + t) - 3.0
    D[2,2]  =  0.0
    D[2,3]  =  4.0*s - 1.0
    D[2,4]  =  0.0
    D[2,5]  = -4.0*r
    D[2,6]  =  4.0*r
    D[2,7]  =  4.0 - 4.0*r - 8.0*s - 4.0*t
    D[2,8]  = -4.0*t
    D[2,9]  =  0.0
    D[2,10] =  4.0*t

    # t-derivatives: dN0/dt to dN9/dt
    D[3,1]  =  4.0*(r + s + t) - 3.0
    D[3,2]  =  0.0
    D[3,3]  =  0.0
    D[3,4]  =  4.0*t - 1.0
    D[3,5]  = -4.0*r
    D[3,6]  =  0.0
    D[3,7]  = -4.0*s
    D[3,8]  =  4.0 - 4.0*r - 4.0*s - 8.0*t
    D[3,9]  =  4.0*r
    D[3,10] =  4.0*s

    return D
end

function shape_func(::Type{Val{WED6}}, R::Array{Float64,1})
    r, s, t = R[1:3]
    N = Array(Float64,6)
    N[1] = 0.5*(1.0-r-s-t+r*t+s*t)
    N[2] = 0.5*(r-r*t)
    N[3] = 0.5*(s-s*t)
    N[4] = 0.5*(1.0-r-s+t-r*t-s*t)
    N[5] = 0.5*(r+r*t)
    N[6] = 0.5*(s+s*t)
    return N
end

function deriv_func(::Type{Val{WED6}}, R::Array{Float64,1})
    r, s, t = R
    D = Array(Float64, 3, 6)
    D[1,1] = 0.5*(-1.0+t);  D[2,1] = 0.5*(-1.0+t);  D[3,1] = 0.5*(-1.0+r+s)
    D[1,2] = 0.5*(1.0-t) ;  D[2,2] = 0.0         ;  D[3,2] = 0.5*(-r)
    D[1,3] = 0.0         ;  D[2,3] = 0.5*(1.0-t) ;  D[3,3] = 0.5*(-s)
    D[1,4] = 0.5*(-1.0-t);  D[2,4] = 0.5*(-1.0-t);  D[3,4] = 0.5*(1.0-r-s)
    D[1,5] = 0.5*(1.0+t) ;  D[2,5] = 0.0         ;  D[3,5] = 0.5*(r)
    D[1,6] = 0.0         ;  D[2,6] = 0.5*(1.0+t) ;  D[3,6] = 0.5*(s)

    return D
end

function shape_func(::Type{Val{WED15}}, R::Array{Float64,1})
    r, s, t = R[1:3]
    N = Array(Float64,15)
    N[1]  = 0.5*(-2.0*(r^2.0)*t-4.0*r*s*t-r*(t^2+0)-2.0*(s^2.0)*t-s*(t^2.0)+2.0*(r^2.0)+4.0*r*s+3.0*r*t+2.0*(s^2.0)+3.0*s*t+(t^2.0)-2.0*r-2.0*s-t) 
    N[2]  = 0.5*(-2.0*(r^2.0)*t+r*(t^2.0)+2.0*(r^2.0)+r*t-2*r)
    N[3]  = 0.5*(-2.0*(s^2.0)*t+s*(t^2.0)+2.0*(s^2.0)+s*t-2.0*s)
    N[4]  = 0.5*(2.0*(r^2.0)*t+4.0*r*s*t-r*(t^2.0)+2.0*(s^2.0)*t-s*(t^2.0)+2.0*(r^2.0)+4.0*r*s-3.0*r*t+2.0*(s^2.0)-3.0*s*t+(t^2.0)-2.0*r-2.0*s+t)
    N[5]  = 0.5*(2.0*(r^2.0)*t+r*(t^2.0)+2.0*(r^2.0)-r*t-2.0*r)
    N[6]  = 0.5*(2.0*(s^2.0)*t+s*(t^2.0)+2.0*(s^2.0)-s*t-2.0*s)
    N[7]  = (2.0*(r^2.0)*t+2*r*s*t-2.0*(r^2.0)-2.0*r*s-2.0*r*t+2*r)
    N[8]  = (-2.0*r*s*t+2.0*r*s)
    N[9]  = (2.0*r*s*t+2.0*(s^2.0)*t-2.0*r*s-2.0*(s^2.0)-2.0*s*t+2.0*s)
    N[10] = (-2.0*(r^2.0)*t-2.0*r*s*t-2.0*(r^2.0)-2.0*r*s+2.0*r*t+2.0*r)
    N[11] = (2.0*r*s*t+2.0*r*s)
    N[12] = (-2.0*r*s*t-2.0*(s^2.0)*t-2.0*r*s-2*(s^2.0)+2.0*s*t+2.0*s)
    N[13] = (r*(t^2.0)+s*(t^2.0)-(t^2.0)-r-s+1.0)
    N[14] = (-r*(t^2.0)+r)
    N[15] = (-s*(t^2.0)+s)
    return N
end

function deriv_func(::Type{Val{WED15}}, R::Array{Float64,1})
    r, s, t = R
    D = Array(Float64, 3, 15)
    # Derivatives with respect to r
    D[1, 1] = 0.5*(-4.0*r*t-4.0*s*t-(t^2.0)+4.0*r+4.0*s+3.0*t-2.0)
    D[1, 2] = 0.5*(-4.0*r*t+(t^2.0)+4.0*r+t-2)
    D[1, 3] = 0.0
    D[1, 4] = 0.5*(4.0*r*t+4.0*s*t-(t^2.0)+4.0*r+4.0*s-3.0*t-2.0)
    D[1, 5] = 0.5*(4.0*r*t+(t^2.0)+4.0*r-t-2.0)
    D[1, 6] = 0.0 
    D[1, 7] = (4.0*r*t+2.0*s*t-4.0*r-2.0*s-2.0*t+2.0) 
    D[1, 8] = (-2.0*s*t+2.0*s)
    D[1, 9] = (2.0*s*t-2.0*s)
    D[1,10] = (-4.0*r*t-2.0*s*t-4.0*r-2.0*s+2.0*t+2.0)
    D[1,11] = (2.0*s*t+2.0*s)
    D[1,12] = (-2.0*s*t-2.0*s)
    D[1,13] = ((t^2.0)-1.0)
    D[1,14] = ((-t^2.0)+1.0)
    D[1,15] = 0.0
   
    # Derivatives with respect to s
    D[2, 1] = 0.5*(-4.0*r*t-4.0*s*t-(t^2.0)+4.0*r+4.0*s+3.0*t-2)
    D[2, 2] = 0.0
    D[2, 3] = 0.5*(-4.0*s*t+(t^2.0)+4*s+t-2.0)
    D[2, 4] = 0.5*(4.0*r*t+4.0*s*t-(t^2.0)+4.0*r+4.0*s-3.0*t-2.0)
    D[2, 5] = 0.0
    D[2, 6] = 0.5*(4.0*s*t+(t^2.0)+4.0*s-t-2.0)
    D[2, 7] = (2.0*r*t-2.0*r)
    D[2, 8] = (-2.0*r*t+2.0*r)
    D[2, 9] = (2.0*r*t+4.0*s*t-2.0*r-4.0*s-2.0*t+2.0)
    D[2,10] = (-2.0*r*t-2.0*r)
    D[2,11] = (2.0*r*t+2.0*r)
    D[2,12] = (-2.0*r*t-4.0*s*t-2.0*r-4.0*s+2.0*t+2.0)
    D[2,13] = ((t^2.0)-1.0)
    D[2,14] = 0.0
    D[2,15] = ((-t^2.0)+1.0) 


    # Derivatives with respect to t
    D[3, 1] = 0.5*(-2.0*(r^2.0)-4.0*r*s-2.0*r*t-2.0*(s^2.0)-2.0*s*t+3.0*r+3.0*s+2.0*t-1.0)
    D[3, 2] = 0.5*(-2.0*(r^2.0)+2.0*r*t+r)
    D[3, 3] = 0.5*(-2.0*(s^2.0)+2*s*t+s)
    D[3, 4] = 0.5*(2.0*(r^2.0)+4.0*r*s-2.0*r*t+2.0*(s^2.0)-2.0*s*t-3.0*r-3.0*s+2.0*t+1)
    D[3, 5] = 0.5*(2.0*(r^2)+2.0*r*t-r)
    D[3, 6] = 0.5*(2*(s^2.0)+2*s*t-s)
    D[3, 7] = (2.0*(r^2.0)+2.0*r*s-2*r)
    D[3, 8] = (-2.0*r*s)
    D[3, 9] = (2.0*r*s+2.0*(s^2.0)-2.0*s)
    D[3,10] = (2.0*r-2.0*(r^2.0)-2.0*r*s)
    D[3,11] = (2.0*r*s)
    D[3,12] = (2.0*s-2.0*s*r-2.0*(s^2.0))
    D[3,13] = (-2.0*t+2.0*t*r+2.0*t*s)
    D[3,14] = (-2.0*r*t)
    D[3,15] = (-2.0*s*t)

    return D
end


# Number of integration points per element
IP_FEM = Dict(
    LIN2    => Dict( 0 => LIN_IP2,  1 => ALL_IP1,  2 => LIN_IP2,  3 => LIN_IP3,   4 => LIN_IP4    ),
    LIN3    => Dict( 0 => LIN_IP2,  2 => LIN_IP2,  3 => LIN_IP3,   4 => LIN_IP4                   ),
    LIN4    => Dict( 0 => LIN_IP3,  2 => LIN_IP2,  3 => LIN_IP3,   4 => LIN_IP4                   ),
    TRI3    => Dict( 0 => TRI_IP1,  1 => TRI_IP1,  3 => TRI_IP3,   6 => TRI_IP6                   ),
    TRI6    => Dict( 0 => TRI_IP3,  3 => TRI_IP3,  6 => TRI_IP6                                   ),
    TRI9    => Dict( 0 => TRI_IP6,  3 => TRI_IP3,  6 => TRI_IP6                                   ),
    TRI10   => Dict( 0 => TRI_IP6,  3 => TRI_IP3,  6 => TRI_IP6                                   ),
    LINK2   => Dict( 0 => LIN_IP2,  2 => LIN_IP2,  3 => LIN_IP3,   4 => LIN_IP4                   ),
    LINK3   => Dict( 0 => LIN_IP3,  2 => LIN_IP2,  3 => LIN_IP3,   4 => LIN_IP4                   ),
    QUAD4   => Dict( 0 => QUAD_IP2, 1 => QUAD_IP1, 4 => QUAD_IP2, 9 => QUAD_IP3, 16 => QUAD_IP4   ),
    QUAD8   => Dict( 0 => QUAD_IP3, 4 => QUAD_IP2, 9 => QUAD_IP3, 16 => QUAD_IP4                  ),
    QUAD9   => Dict( 0 => QUAD_IP3, 4 => QUAD_IP2, 9 => QUAD_IP3, 16 => QUAD_IP4                  ),
    QUAD12  => Dict( 0 => QUAD_IP4, 4 => QUAD_IP2, 9 => QUAD_IP3, 16 => QUAD_IP4                  ),
    QUAD16  => Dict( 0 => QUAD_IP4, 4 => QUAD_IP2, 9 => QUAD_IP3, 16 => QUAD_IP4                  ),
    TET4    => Dict( 0 => TET_IP4,  1 => TET_IP1,  4 => TET_IP4,   5 => TET_IP5,  11 => TET_IP11  ),
    TET10   => Dict( 0 => TET_IP4,  1 => TET_IP1,  4 => TET_IP4,   5 => TET_IP5,  11 => TET_IP11  ),
    HEX8    => Dict( 0 => HEX_IP2,  8 => HEX_IP2, 27 => HEX_IP3                                   ),
    HEX20   => Dict( 0 => HEX_IP3,  8 => HEX_IP2, 27 => HEX_IP3                                   ),
    WED6    => Dict( 0 => WED_IP9,  2 => WED_IP2,  9 => WED_IP9,  18 => WED_IP18                  ),
    WED15   => Dict( 0 => WED_IP9,  2 => WED_IP2,  9 => WED_IP9,  18 => WED_IP18                  ),

    JLIN2   => Dict( 0 => LIN_IP2,  2 => LIN_IP2,  3 => LIN_IP3,   4 => LIN_IP4                   ),
    JLIN3   => Dict( 0 => LIN_IP2,  2 => LIN_IP2,  3 => LIN_IP3,   4 => LIN_IP4                   ),
    JLIN4   => Dict( 0 => LIN_IP3,  2 => LIN_IP2,  3 => LIN_IP3,   4 => LIN_IP4                   ),
    JTRI3   => Dict( 0 => TRI_IP3,  3 => TRI_IP3,  6 => TRI_IP6                                   ),
    JTRI6   => Dict( 0 => TRI_IP3,  3 => TRI_IP3,  6 => TRI_IP6                                   ),
    JQUAD4  => Dict( 0 => QUAD_IP2, 4 => QUAD_IP2, 9 => QUAD_IP3, 16 => QUAD_IP4                  ),
    JQUAD8  => Dict( 0 => QUAD_IP3, 4 => QUAD_IP2, 9 => QUAD_IP3, 16 => QUAD_IP4                  )
)

function get_ip_coords(shape::ShapeType, nips=0)
    if !haskey(IP_FEM, shape)
        error("Ip coordinates for shape ($shape) is not available")
    end
    all_shape_coord = IP_FEM[shape]
    if !haskey(all_shape_coord, nips)
        error("Number of ips ($nips) for shape ($shape) is not available")
    end
    all_shape_coord[nips]
end

function bdistance(shape::ShapeType, R::Array{Float64,1})
    # Returns a real value which is a pseudo distance from a point to the border of an element
    # Arguments:
    #     R - a vector containing the point coordinates
    # Returns:
    #     a real value: if possitive then the point is inside the element and negative otherwise
    
    r, s, t = R
    if shape in [ TRI3, TRI6, TRI9, TRI10 ] return min(r, s, 1.0-r-s) end
    if shape in [ QUAD4, QUAD8, QUAD12, QUAD16 ]  return min(1.0 - abs(r), 1.0 - abs(s)) end
    if shape in [ TET4, TET10 ] return min(r, s, t, 1.0-r-s-t) end
    if shape in [ HEX8, HEX20 ] return min(1.0 - abs(r), 1.0 - abs(s), 1.0 - abs(t)) end
    if shape in [ WED6, WED15 ] return min(r, s, 1.0-r-s, 1.0-abs(t)) end
    #if shape in [ QUAD4, QUAD8, QUAD12, QUAD16 ]  return min(1.0 - r*r, 1.0 - s*s) end
    #if shape in [ HEX8, HEX20 ] return min(1.0 - r*r, 1.0 - s*s, 1.0 - t*t) end
    error("No boundary distance for shape ($shape)")
end

function inverse_map(shape::ShapeType, coords::Array{Float64,2}, X0::Array{Float64,1}, Tol=1.0e-7)
    MAXIT = 20
    ndim  = get_ndim(shape)
    R = zeros(ndim)
    C = coords
    local ΔX::Array{Float64,1}

    X = X0
    if size(coords,2)==2
        X = X0[1:ndim]
    end

    k = 0
    for k=1:MAXIT
        # calculate Jacobian
        D = deriv_func(shape, R)
        J = D*C

        # calculate trial of real coordinates
        N  = shape_func(shape, R)
        Xt = C'*N # interpolating

        # calculate the error
        ΔX = Xt - X
        ΔR = pinv(J)'*ΔX

        # updating local coords R
        R -= ΔR
        if norm(ΔX) < Tol; break end
    end

    # TODO: Improve accuracy of inverse_map function in elements with non regular shape
    #k==MAXIT && println("Warning: max iterations (MAXIT=$MAXIT) reached in inverse mapping. norm(ΔX)=$(norm(ΔX))")

    if ndim==2
        R = vcat( R, 0.0 )
    end
    return R
end


function is_inside(shape::ShapeType, C::Array{Float64,2}, X::Array{Float64,1}, Tol = 1.e-7)
    if !is_solid(shape) return false end

    # Testing with bounding box
    ndim = size(C,1)
    Cmin = vec(minimum(C,1))
    Cmax = vec(maximum(C,1))
    maxl = maximum(Cmax-Cmin)
    ttol = 0.1*maxl # 10% is important for curved elements

    if any(X .< Cmin-ttol) || any(X .> Cmax+ttol)
        return false
    end

    # Testing with inverse mapping
    R = inverse_map(shape, C, X, Tol)
    if bdistance(shape, R) > -Tol
        return true
    else
        return false
    end
end

function extrapolator(shape::ShapeType, nips::Int)
    #  Returns a numpy matrix E that extrapolates ip values to nodal values as:
    #
    #                 NodalValues = E * IpValues;
    # where:
    #                             +                +              +
    #                        E = N * (I - EPS1*EPS1 ) + EPS * EPS1
    #
    # and            N = [shape functions matrix]
    #                         1        2        ...        nNodes
    #                  1    [[N_11 N_12
    #                  2     [N_21
    #                  :     [
    #                 nIP    [N_ ...                    ]]
    

    nnodes  = get_nnodes(shape)
    IP      = get_ip_coords(shape, nips)
    ndim    = get_ndim(shape) # shape ndim: not related with the analysis ndim

    #filling N matrix with shape functions of all ips
    N = Array(Float64, nips, nnodes)
    for i=1:nips
        N[i,:] = shape_func(shape, vec(IP[i,:]))
    end

    #calculate extrapolator matrix
    if nips==nnodes
        return inv(N)
    elseif nips>=nnodes
        return pinv(N)
    elseif nips==1
        return pinv(N) # Correction procedure is not applicable for nips==1
    end

    I = eye(nips)

    # εip matrix: Local ip coordinates of integration points
    εip = [ IP[:,1:ndim] zeros(nips) ]

    #@show nips
    #@show nnodes
    # ε matrix: Local coordinates of nodal points
    ε = get_local_coords(shape)

    E = pinv(N)*(I - εip*pinv(εip)) + ε*pinv(εip)

    return E
end
