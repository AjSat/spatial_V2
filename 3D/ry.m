function  E = ry( theta )

% ry  3x3 coordinate rotation (Y-axis)
% ry(theta)  calculates the 3x3 rotational coordinate transform matrix from
% A to B coordinates, where coordinate frame B is rotated by an angle theta
% (radians) relative to frame A about their common Y axis.

import casadi.*
c = cos(theta);
s = sin(theta);

if strcmp(class(theta), 'casadi.MX')
    E = MX(3,3);
else
    E = SX(3,3);
end
E(1,1) = c;
E(2,2) = 1;
E(3,3) = c;
E(1,3) = -s;
E(3,1) = s;

% E = [ c  0 -s;
%       0  1  0;
%       s  0  c ];
