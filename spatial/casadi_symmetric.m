function  M = casadi_symmetric( M )

% Takes a symmetric casadi matrix and duplicates the symbols.

m = size(M,1);

for i = 2:m
    for j = 1:i-1
        
        M(i,j) = M(j,i);
    end
end

