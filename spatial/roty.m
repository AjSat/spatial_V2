function  X = roty( theta )

% roty  spatial coordinate transform (Y-axis rotation).
% roty(theta)  calculates the coordinate transform matrix from A to B
% coordinates for spatial motion vectors, where coordinate frame B is
% rotated by an angle theta (radians) relative to frame A about their
% common Y axis.

import casadi.*
roty3 = ry(theta);

if strcmp(class(theta), 'casadi.MX')
    X = MX(6,6);
else
    X = SX(6,6);
end
X(1:3, 1:3) = roty3;
X(4:6, 4:6) = roty3;

% X = [ c  0 -s  0  0  0 ;
%       0  1  0  0  0  0 ;
%       s  0  c  0  0  0 ;
%       0  0  0  c  0 -s ;
%       0  0  0  0  1  0 ;
%       0  0  0  s  0  c
%     ];
