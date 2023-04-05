function [X_inv] = spatial_inv(X)
%SPATIAL_INV Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
X_inv = SX(6,6);
X_inv(1:3, 1:3) = X(1:3, 1:3)';
X_inv(4:6, 4:6) = X(1:3, 1:3)';
X_inv(4:6, 1:3) = -X(1:3, 1:3)'*X(4:6, 1:3)*X(1:3, 1:3)';
end

