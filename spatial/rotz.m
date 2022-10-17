function  X = rotz( theta )

% rotz  spatial coordinate transform (Z-axis rotation).
% rotz(theta)  calculates the coordinate transform matrix from A to B
% coordinates for spatial motion vectors, where coordinate frame B is
% rotated by an angle theta (radians) relative to frame A about their
% common Z axis.

import casadi.*
rotz3 = rz(theta);

if strcmp(class(theta), 'casadi.MX')
    X = MX(6,6);
else
    X = SX(6,6);
end
X(1:3, 1:3) = rotz3;
X(4:6, 4:6) = rotz3;

% X = [  c  s  0  0  0  0 ;
%       -s  c  0  0  0  0 ;
%        0  0  1  0  0  0 ;
%        0  0  0  c  s  0 ;
%        0  0  0 -s  c  0 ;
%        0  0  0  0  0  1
%     ];
