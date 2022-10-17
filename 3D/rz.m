function  E = rz( theta )

% rz  3x3 coordinate rotation (Z-axis)
% rz(theta)  calculates the 3x3 rotational coordinate transform matrix from
% A to B coordinates, where coordinate frame B is rotated by an angle theta
% (radians) relative to frame A about their common Z axis.

import casadi.*
c = cos(theta);
s = sin(theta);

if strcmp(class(theta), 'casadi.MX')
    E = MX(3,3);
else
    E = SX(3,3);
end
E(1,1) = c;
E(2,2) = c;
E(3,3) = 1;
E(1,2) = s;
E(2,1) = -s;
% E = [ c  s  0;
%      -s  c  0;
%       0  0  1 ];
