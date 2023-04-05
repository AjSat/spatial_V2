function [X12] = compose_frames(X1,X2)
%COMPOSE_FRAMES Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
X12 = SX(6,6);
X12(1:3,1:3) = X1(1:3,1:3)*X2(1:3,1:3);
X12(4:6,4:6) = X12(1:3,1:3);
X12(4:6, 1:3) = X1(4:6,1:3)*X2(1:3,1:3) + X1(4:6,4:6)*X2(4:6,1:3);

end

