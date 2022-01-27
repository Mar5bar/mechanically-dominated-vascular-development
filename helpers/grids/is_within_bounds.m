function r = is_within_bounds(x,y,params)
% IS_WITHIN_BOUNDS checks if the point (x,y) is inside of the domain.
% Note that the domain is completely specified in this function.
    if nargin == 1
        y = x(:,2);
        x = x(:,1);
    end

    % Annulus, radii 0.1 and 1.
    rad = sqrt(x.^2 + y.^2);
    r = (rad < params.outerRadius) & (rad > params.innerRadius);
end