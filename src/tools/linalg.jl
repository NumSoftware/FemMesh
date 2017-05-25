
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

# Types
const Vect = Array{Float64, 1}
const Matx = Array{Float64, 2}

# Fancy matrix printing
function print_matrix(M::Array{Float64,2})
    n, m = size(M)
    for i=1:n
        for j=1:m
            @printf( "%23.11e", M[i,j] )
        end
        println()
    end
end

# Pseudo determinant of non-square matrices
function norm2(J)

    if ndims(J)==1; return norm(J) end

    r, c = size(J)
    if r==1; return norm(J) end
    if r==2 && c==3
        j1 = J[1,1]*J[2,2] - J[1,2]*J[2,1]
        j2 = J[1,2]*J[2,3] - J[1,3]*J[2,2]
        j3 = J[1,3]*J[2,1] - J[1,1]*J[2,3]
        return (j1*j1 + j2*j2 + j3*j3)^0.5  # jacobian determinant
    end
    if r==c; return det(J) end
    error("No rule to calculate norm2 of a $r x $c matrix")
end


macro gemm(expr)
    β = 0.0
    s = 1.0
    if expr.head == :(=)
    elseif expr.head == :(+=)
        β = 1.0
    elseif expr.head == :(-=)
        s = -1.0
        β =  1.0
    else
        error("@gemm: =, +=, -= operator expected, found $(expr.head)")
    end

    C = expr.args[1]
    rhs = expr.args[2]

    if rhs.args[1] !=  :(*)
        error("@inplace: * operator expected, found $(rhs.args[1])")
    end

    if length(rhs.args) == 4
        α = rhs.args[2]
        A = rhs.args[3]
        B = rhs.args[4]
    else
        α = 1.0
        A = rhs.args[2]
        B = rhs.args[3]
    end
    
    tA = 'N'
    if typeof(A) == Expr 
        if A.head == Symbol("'"); 
            tA = 'T' 
            A  = A.args[1]
        end
    end

    tB = 'N'
    if typeof(B) == Expr
        if B.head == Symbol("'"); 
            tB = 'T' 
            B  = B.args[1]
        end
    end

    return :( BLAS.gemm!($tA, $tB, $(esc(α))*$s, $(esc(A)), $(esc(B)), $β, $(esc(C)) ) )
end


macro gemv(expr)
    β = 0.0
    s = 1.0
    if expr.head == :(=)
    elseif expr.head == :(+=)
        β = 1.0
    elseif expr.head == :(-=)
        s = -1.0
        β =  1.0
    else
        error("@inplace: =, +=, -= operator expected, found $(expr.head)")
    end

    C = expr.args[1]
    rhs = expr.args[2]

    if rhs.args[1] !=  :(*)
        error("@inplace: * operator expected, found $(rhs.args[1])")
    end

    if length(rhs.args) == 4
        α = rhs.args[2]
        A = rhs.args[3]
        B = rhs.args[4]
    else
        α = 1.0
        A = rhs.args[2]
        B = rhs.args[3]
    end
    
    tA = 'N'
    if typeof(A) == Expr 
        if A.head == Symbol("'"); 
            tA = 'T' 
            A  = A.args[1]
        end
    end

    return :( BLAS.gemv!($tA, $(esc(α))*$s, $(esc(A)), $(esc(B)), $(β), $(esc(C)) ) )
end

# Y += α*X   
# Y  = α*X
macro scale(expr)
    α = 1.0
    β = 0.0
    s = 1.0
    if expr.head == :(=)
    elseif expr.head == :(+=)
        β = 1.0
    elseif expr.head == :(-=)
        s = -1.0
        β =  1.0
    else
        error("@scale: =, +=, -= operator expected, found $(expr.head)")
    end

    Y   = expr.args[1]
    rhs = expr.args[2]

    if typeof(rhs)==Expr
        if rhs.args[1] !=  :(*)
            error("@scale: * operator expected, found $(rhs.args[1])")
        end

        if length(rhs.args) == 3
            α = rhs.args[2]
            X = rhs.args[3]
        else
            X = rhs.args[2]
        end
    else
        X = rhs
    end

    if β == 0.0
        return quote
            $(esc(Y))[:] = 0.0
            BLAS.axpy!( $(esc(α))*$s, $(esc(X)), $(esc(Y)) )
        end
    else
        return :( BLAS.axpy!( $(esc(α))*$s, $(esc(X)), $(esc(Y)) ) )
    end
end


