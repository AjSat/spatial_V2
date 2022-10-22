function [x] = back_sub(U,v)

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

[nRow,nCol]=size(v);
x=SX(nRow,nCol);


for c = 1:nCol
    x(nRow, c)=v(nRow, c)/U(nRow,nRow);
    for n=(nRow-1):-1:1
        x(n, c)=(v(n, c)-(U(n,n+1:end)*x(n+1:end, c))) / U(n,n);
    end
end