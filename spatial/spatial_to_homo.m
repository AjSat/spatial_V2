function [o1] = spatial_to_homo(X)
%SPATIAL_TO_HOMO Summary of this function goes here
%   Detailed explanation goes here
o1 = eye(4,4);
o1(1:3,1:3) = X(1:3,1:3)';
vec_skew = X(1:3,1:3)'*X(4:6,1:3);
o1(1,4) = vec_skew(2, 3);
o1(2,4) = vec_skew(3, 1);
o1(3,4) = vec_skew(1,2);
end

