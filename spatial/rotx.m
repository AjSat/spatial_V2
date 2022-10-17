function  X = rotx( theta )

% rotx  spatial coordinate transform (X-axis rotation).
% rotx(theta)  calculates the coordinate transform matrix from A to B
% coordinates for spatial motion vectors, where coordinate frame B is
% rotated by an angle theta (radians) relative to frame A about their
% common X axis.

import casadi.*
rotx3 = rx(theta);

if strcmp(class(theta), 'casadi.MX')
    X = MX(6,6);
else
    X = SX(6,6);
end

X(1:3, 1:3) = rotx3;
X(4:6, 4:6) = rotx3;

% X = [ 1  0  0  0  0  0 ;
%       0  c  s  0  0  0 ;
%       0 -s  c  0  0  0 ;
%       0  0  0  1  0  0 ;
%       0  0  0  0  c  s ;
%       0  0  0  0 -s  c
%     ];
