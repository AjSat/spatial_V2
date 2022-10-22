function [y] = forward_sub(L,w)

%Function to find the Cholesky factor of a Positive Definite matrix A
%Author Mathuranathan for https://www.gaussianwaves.com
%Licensed under Creative Commons: CC-NC-BY-SA 3.0
%A = positive definite matrix
%Option can be one of the following 'Lower','Upper'
%L = Cholesky factorizaton of A such that A=LL^T
%If option ='Lower', then it returns the Cholesky factored matrix L in
%lower triangular form
%If option ='Upper', then it returns the Cholesky factored matrix L in
%upper triangular form

%Test for positive definiteness (symmetricity need to satisfy)
%Check if the matrix is symmetric

import casadi.*
reg = 1e-12;

[nRow,nCol]=size(w);
y= SX(nRow,nCol);


for c = 1:nCol
    y(1, c)=w(1, c)/L(1,1);
    for n=2:nRow
        y(n, c)=( w(n, c) - L(n,1:n-1)*y(1:n-1, c) ) / L(n,n);
    end
end