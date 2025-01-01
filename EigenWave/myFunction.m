function [U] = myFunction(W0,W1,par,fun)
W0 = Vector_to_Matrix(W0,par);
W0 = expand_matrix(W0,par);
U = WaveHoltz(W0,W1,par,fun);
U = shrink_matrix(U,par);
U = Matrix_to_Vector(U);
end