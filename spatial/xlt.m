function  X = xlt( r )

% xlt  spatial coordinate transform (translation of origin).
% xlt(r)  calculates the coordinate transform matrix from A to B
% coordinates for spatial motion vectors, in which frame B is translated by
% an amount r (3D vector) relative to frame A.  r can be a row or column
% vector.

% X = [  1     0     0    0  0  0 ;
%        0     1     0    0  0  0 ;
%        0     0     1    0  0  0 ;
%        0     r(3) -r(2) 1  0  0 ;
%       -r(3)  0     r(1) 0  1  0 ;
%        r(2) -r(1)  0    0  0  1
%     ];

import casadi.*;

if strcmp(class(r), 'casadi.MX')
    X = MX(6,6);
else
    r = sparsify(SX(r));
    X = SX(6,6);
end
for i = 1:6
    X(i,i) = 1;
end

X(4, 2) = r(3);
X(5, 1) = -r(3);
X(4, 3) = -r(2);
X(6, 1) = r(2);
X(6, 2) = -r(1);
X(5, 3) = r(1);



