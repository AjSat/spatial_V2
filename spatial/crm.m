function  vcross = crm( v )

% crm  spatial/planar cross-product operator (motion).
% crm(v)  calculates the 6x6 (or 3x3) matrix such that the expression
% crm(v)*m is the cross product of the motion vectors v and m.  If
% length(v)==6 then it is taken to be a spatial vector, and the return
% value is a 6x6 matrix.  Otherwise, v is taken to be a planar vector, and
% the return value is 3x3.

import casadi.*

if length(v) == 6
    
  if strcmp(class(v), 'casadi.MX')
    vcross = MX(6,6);
  else
    vcross = SX(6,6);
  end
%   vcross = SX(6,6);
  vcross(1, 2) = -v(3);
  vcross(2, 1) = v(3);
  vcross(1, 3) = v(2);
  vcross(3, 1) = -v(2);
  vcross(3, 2) = v(1);
  vcross(2, 3) = -v(1);
  
  vcross(4:6, 4:6) = vcross(1:3, 1:3);
  
  vcross(4, 2) = -v(6);
  vcross(5, 1) = v(6);
  vcross(4, 3) = v(5);
  vcross(6, 1) = -v(5);
  vcross(6, 2) = v(4);
  vcross(5, 3) = -v(4);
%   vcross(4:6, 1:3) = [0    -v(6)  v(5);
% 	      v(6)  0    -v(4);
% 	     -v(5)  v(4)  0 ];
%   vcross = [  0    -v(3)  v(2)   0     0     0    ;
% 	      v(3)  0    -v(1)   0     0     0    ;
% 	     -v(2)  v(1)  0      0     0     0    ;
% 	      0    -v(6)  v(5)   0    -v(3)  v(2) ;
% 	      v(6)  0    -v(4)   v(3)  0    -v(1) ;
% 	     -v(5)  v(4)  0     -v(2)  v(1)  0 ];

else

  vcross = [  0     0     0    ;
	      v(3)  0    -v(1) ;
	     -v(2)  v(1)  0 ];
end
