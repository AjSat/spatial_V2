function [y] = forward_sub(L,w)

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