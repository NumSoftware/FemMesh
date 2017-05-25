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


export ShapeType
export ShapeClass
export get_ip_coords, get_shape_from_vtk
export inverse_map, extrapolator

@enum(ShapeClass,
LINE_SHAPE    = 1,
SOLID_SHAPE   = 2,
JOINT_SHAPE   = 3,
JOINT1D_SHAPE = 4,
EMBEDDED      = 5
)

# Export
for s in instances(ShapeClass)
    @eval export $(Symbol(s))
end 



type ShapeType 
    name       ::String
    class      ::ShapeClass
    ndim       ::Int
    npoints    ::Int
    basic_shape::ShapeType
    vtk_type   ::Int
    facet_idxs ::Array
    edge_idxs  ::Array
    facet_shape::Union{ShapeType, Tuple}
    nat_coords ::Array
    quadrature ::Dict{Int, Array}
    func       ::Function
    deriv      ::Function
    function ShapeType()
        return new()
    end
end

# Shape for unknown polyvertex
const POLYV = ShapeType()


include("cells/quadrature.jl")
include("cells/vtk.jl")
include("cells/lines.jl")
include("cells/solids2d.jl")
include("cells/solids3d.jl")
include("cells/joints.jl")


function get_ip_coords(shape::ShapeType, nips=0)::Array{Float64,2}
    !haskey(shape.quadrature, nips) && error("Number of ips ($nips) for shape ($shape) is not available")
    return shape.quadrature[nips]
end


# Available shape types
const VTK2SHAPE = Dict{VTKCellType,ShapeType}(
    VTK_POLY_VERTEX                      => POLYV,
    VTK_LINE                             => LIN2,
    VTK_TRIANGLE                         => TRI3,
    VTK_QUAD                             => QUAD4,
    VTK_TETRA                            => TET4,
    VTK_HEXAHEDRON                       => HEX8,
    VTK_WEDGE                            => WED6,
    VTK_QUADRATIC_EDGE                   => LIN3,
    VTK_QUADRATIC_TRIANGLE               => TRI6,
    VTK_QUADRATIC_QUAD                   => QUAD8,
    VTK_QUADRATIC_TETRA                  => TET10,
    VTK_QUADRATIC_HEXAHEDRON             => HEX20,
    VTK_QUADRATIC_WEDGE                  => WED15,
    VTK_BIQUADRATIC_QUAD                 => QUAD9
)


function get_shape_from_vtk(vtk_type::VTKCellType, npoints::Int64, ndim::Int64)::ShapeType

    if vtk_type!=VTK_POLY_VERTEX return VTK2SHAPE[vtk_type] end

    if     npoints==2   return JLINK2
    elseif npoints==3   return JLINK3
    elseif npoints==4   return JLIN2
    elseif npoints==6   
        if ndim==2 return JLIN3 end
        if ndim==3 return JTRI3 end
    elseif npoints==8
        if ndim==2 return JLIN4  end
        if ndim==3 return JQUAD4 end
    elseif npoints==9   return TRI9
    elseif npoints==10  return TRI10
    elseif npoints==12  
        #if ndim==2 return QUAD12 end
        if ndim==3 return JTRI6  end
    elseif npoints==16
        if ndim==2 return QUAD16 end
        if ndim==3 return JQUAD8 end
    #elseif npoints==18  return JTRI9
    #elseif npoints==20  return JTRI10
    end

    error("get_shape_from_vtk: Unknown shape for vtk_type $vtk_type and npoints $npoints with ndim $ndim")
end


function bdistance(shape::ShapeType, R::Array{Float64,1})
    # Returns a real value which is a pseudo distance from a point to the border of an element
    # Arguments:
    #     R - a vector containing the point coordinates
    # Returns:
    #     a real value: if possitive then the point is inside the element and negative otherwise
    
    r, s, t = R
    bshape = shape.basic_shape
    if bshape == TRI3  return min(r, s, 1.0-r-s) end
    if bshape == QUAD4 return min(1.0 - abs(r), 1.0 - abs(s)) end
    if bshape == TET4  return min(r, s, t, 1.0-r-s-t) end
    if bshape == HEX8  return min(1.0 - abs(r), 1.0 - abs(s), 1.0 - abs(t)) end
    if bshape == WED6  return min(r, s, 1.0-r-s, 1.0-abs(t)) end
    error("No boundary distance for shape ($shape)")
end


function inverse_map(shape::ShapeType, coords::Array{Float64,2}, X0::Array{Float64,1}, Tol=1.0e-7)
    MAXIT = 20
    ndim  = shape.ndim
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
        D = shape.deriv(R)
        J = D*C

        # calculate trial of real coordinates
        N  = shape.func(R)
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
    if shape.class!=SOLID_SHAPE return false end

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
    

    npoints = shape.npoints
    IP      = get_ip_coords(shape, nips)
    ndim    = shape.ndim # not related with the analysis ndim

    #filling N matrix with shape functions of all ips
    N = Array{Float64}(nips, npoints)
    for i=1:nips
        N[i,:] = shape.func(vec(IP[i,:]))
    end

    #calculate extrapolator matrix
    if nips==npoints
        return inv(N)
    elseif nips>=npoints
        return pinv(N)
    elseif nips==1
        return pinv(N) # Correction procedure is not applicable for nips==1
    end

    I = eye(nips)

    # εip matrix: Local ip coordinates of integration points
    εip = [ IP[:,1:ndim] ones(nips) ]

    #@show nips
    #@show npoints
    # ε matrix: Local coordinates of nodal points
    ε = [ get_local_coords(shape) ones(npoints) ] # increase a column of ones

    E = pinv(N)*(I - εip*pinv(εip)) + ε*pinv(εip)

    return E
end
