function [x] = back_sub(U,v)

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