function [agent_indices, agent_positions] =  init_agent_positions(X_grid,Y_grid,N_agents)
% INIT_WALKER_POSITIONS generates and returns the grid indices and raw positions of N_agents
% initial agents.
    x = zeros(N_agents,1); y = zeros(N_agents,1);
    % Uniformly spaced on a ring of radius 0.1.
    phis = linspace(0,2*pi,N_agents+1)'; phis = phis(1:N_agents);
    for i = 1 : N_agents
        phi = phis(i);
        x(i) = 0.11 * cos(phi);
        y(i) = 0.11 * sin(phi);
    end
    [x,y,i,j] = position_to_grid(x,y,X_grid,Y_grid);
    agent_indices = [i(:),j(:)]; agent_positions = [x(:),y(:)];
end