function basis = basis_func(x,y,params)
%% Evaluate the basis function at all pairs of points given in the arrays.
    basis = 1./(params.lambda*((x-x').^2 + (y-y').^2)+1);
end