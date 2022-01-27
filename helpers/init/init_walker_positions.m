function [walker_indices, walker_positions] =  init_walker_positions(X_grid,Y_grid,N_walkers)
% INIT_WALKER_POSITIONS generates and returns the grid indices and raw positions of N_walkers
% initial walkers.
    x = zeros(N_walkers,1); y = zeros(N_walkers,1);
    % Uniformly spaced on a ring of radius 0.1.
    phis = linspace(0,2*pi,N_walkers+1)'; phis = phis(1:N_walkers);
    for i = 1 : N_walkers
        phi = phis(i);
        x(i) = 0.11 * cos(phi);
        y(i) = 0.11 * sin(phi);
    end
    [x,y,i,j] = position_to_grid(x,y,X_grid,Y_grid);
    walker_indices = [i(:),j(:)]; walker_positions = [x(:),y(:)];
end