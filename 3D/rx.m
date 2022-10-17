function  E = rx( theta )

% rx  3x3 coordinate rotation (X-axis)
% rx(theta)  calculates the 3x3 rotational coordinate transform matrix from
% A to B coordinates, where coordinate frame B is rotated by an angle theta
% (radians) relative to frame A about their common X axis.

import casadi.*
c = cos(theta);
s = sin(theta);

if strcmp(class(theta), 'casadi.MX')
    E = MX(3,3);
else
    E = SX(3,3);
end
E(1,1) = 1;
E(2,2) = c;
E(3,3) = c;
E(2,3) = s;
E(3,2) = -s;
% E = [ 1  0  0;
%       0  c  s;
%       0 -s  c ];
