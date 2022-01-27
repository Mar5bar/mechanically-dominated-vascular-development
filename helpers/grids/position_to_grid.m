function [x,y,i,j,e] = position_to_grid(x_raw,y_raw,X_grid,Y_grid)
% POSITION_TO_GRID returns the positions x_raw, y_raw rounded on to the grid with
% x coords X_grid and y coords Y_grid (these are vectors). Accepts vector x,y.
    [ex,i] = min(abs(x_raw(:) - X_grid),[],2);
    x = X_grid(i);
    [ey,j] = min(abs(y_raw(:) - Y_grid),[],2);
    y = Y_grid(j);
    if nargout > 4
        e = (ex.^2 + ey.^2).^0.5;
    end
end