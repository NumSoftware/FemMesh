include("../tools/linalg.jl")

A=[1. 2 3; 4 5 6; 1 2 3 ]
B=[1. 1 1; 2 2 2; 3 3 3 ]
A = [ 1., 1, 1 ]
B = [ 1., 2, 3 ]
C=ones(3)

#@gemv C += 2*A'*B

@scale A += 2*B

#@show C
#@show ones(3) + 2*A'*B - C
@show ones(3) + 2*B - A

#@show macroexpand( (:( @gemv C-=coef*A'*B  ) ))
@show macroexpand( (:( @scale A-=2*B ) ))
