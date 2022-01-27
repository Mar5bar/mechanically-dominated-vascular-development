function [i,j] = coordinate_to_index(x,y,X_grid,Y_grid)
    i = find(X_grid==x);
    j = find(Y_grid==y);
end