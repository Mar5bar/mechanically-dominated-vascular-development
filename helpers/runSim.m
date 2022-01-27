function [trails, phi, grids, event_log, basis] = runSim(params)

    % Initialise.
    grid_X = linspace(params.grid_X_lims(1),params.grid_X_lims(2),params.grid_X_num);
    grid_Y = linspace(params.grid_Y_lims(1),params.grid_Y_lims(2),params.grid_Y_num);
    [mesh_X,mesh_Y] = meshgrid(grid_X,grid_Y); % These grids are transposed to row-column indexing.
    mesh_X = transpose(mesh_X); mesh_Y = transpose(mesh_Y);
    trails = sparse([],[],[],params.grid_X_num,params.grid_Y_num,50000);
    % Generate a (potentially distinct) grid for phi.
    phi_grid_X = linspace(params.grid_X_lims(1),params.grid_X_lims(2),params.phi_grid_X_num);
    phi_grid_Y = linspace(params.grid_Y_lims(1),params.grid_Y_lims(2),params.phi_grid_Y_num);
    [phi_mesh_X,phi_mesh_Y] = meshgrid(phi_grid_X,phi_grid_Y); % These grids are transposed to row-column indexing.
    phi_mesh_X = transpose(phi_mesh_X); phi_mesh_Y = transpose(phi_mesh_Y);
    within_bounds_phi_mask = sparse(is_within_bounds(phi_mesh_X,phi_mesh_Y,params));
    numbered_phi_grid = double(within_bounds_phi_mask);
    numbered_phi_grid(within_bounds_phi_mask) = 1:length(find(numbered_phi_grid));
    numRand = params.N_agents*params.N_steps; % Number of random numbers/choices to generate at once.

    % Structures for saving solutions after each repeat.
    weights_each_rep = cell(1,params.num_repeats);
    trails_each_rep = cell(1,params.num_repeats);
    event_log = cell(1,params.num_repeats);

    % Precompute the basis function evaluated on the mesh.
    x = phi_mesh_X(:); y = phi_mesh_Y(:);
    basis = basis_func(x,y,params);

    if params.num_repeats > 1
        textprogressbar('Progress: ',0)
    end
    
    for rep_num = 1 : params.num_repeats
        N_agents = params.N_agents;
        trails(:) = 0;
        [agent_indices, agent_positions] = init_agent_positions(grid_X,grid_Y,N_agents);
        trails = draw_agents(1:N_agents,agent_indices,trails);
        is_active = true(1,N_agents);
        agent_directions = init_agent_directions(agent_positions);
        within_bounds_mask = sparse(is_within_bounds(mesh_X,mesh_Y,params));
        splitting_counters = init_splitting_counters(params.splitting_period,N_agents); % Steps until splitting.

        %-----------
        % For fast evaluation of phi:
        if params.reset_phi || (~exist('weights'))
            phi_0 = init_phi(phi_mesh_X,phi_mesh_Y,within_bounds_phi_mask,params.init_phi);
            weights = zeros(1,params.phi_grid_X_num*params.phi_grid_Y_num);
        elseif params.degrade_phi
            % Reduce the magnitude of the modifications to the field by 1/2.
            weights = weights / 2;
        end
        % We can now evaluate via:
        % sum_T = dot(weights,basis(:,ind))';
        % where ind is the linear index of the evaluation point.
        % Evaluation at multiple points is done via 
        % sum_T = weights * basis(:,inds)';

        %-----------

        % Generate random numbers. For collisions/merge events, we will need a
        % large number of random samples from {-1,1} and random numbers.
        [rand_choices, rand_direction, rand_choice_counter, rand_direction_counter] = generateRandoms(numRand,params);

        % Arrays for logging events.
        merges = [];
        splits = [];
        collisions = [];

        % Step
        for step_ind = 1 : params.N_steps

            % Random order of agents.
            active_agents = shuffle(find(is_active));

            % Evaluate phi at the locations of the active agents.
            % Round the agents to the phi grid and use the nearest gridpoint.
            [x,y,i,j] = position_to_grid(agent_positions(active_agents,1),agent_positions(active_agents,2),phi_grid_X,phi_grid_Y);
            lin_inds = i + (j-1) * params.phi_grid_X_num;
            phi_at_agents = eval_phi_basis(lin_inds, phi_0(lin_inds(:)), basis, weights);
            phi_at_agents(~within_bounds_phi_mask(lin_inds)) = NaN;

            % If we don't have enough random numbers, generate more.
            if rand_direction_counter + length(active_agents) > numRand
                [rand_choices, rand_direction, rand_choice_counter, rand_direction_counter] = generateRandoms(numRand,params);
            end
            % Move agents, recording trails and positions.
            for ind = 1:length(active_agents)
                w = active_agents(ind);
                agent_directions(w) = gen_agent_direction(w,agent_directions(w),phi_at_agents(ind),rand_direction(rand_direction_counter),params.k);
                rand_direction_counter = rand_direction_counter + 1;
                % Update the (unrounded) position.
                agent_positions(w,:) = agent_positions(w,:) + params.speed * [cos(agent_directions(w)),sin(agent_directions(w))];
            end
            [x,y,i,j] = position_to_grid(agent_positions(active_agents,1),agent_positions(active_agents,2),grid_X,grid_Y);
            % Draw the paths and check collisions for all agents.
            [trails,new_collisions,endpoints] = draw_paths(active_agents,agent_indices(active_agents,:),[i,j],trails);

            % Add new weights.
            % Convert the start and endpoints of the paths (indices) to phi_indices.
            [~,~,istart,jstart] = position_to_grid(grid_X(agent_indices(active_agents,1)),grid_Y(agent_indices(active_agents,2)),phi_grid_X,phi_grid_Y);
            [~,~,iend,jend] = position_to_grid(grid_X(endpoints(:,1)),grid_Y(endpoints(:,2)),phi_grid_X,phi_grid_Y);
            % For each agent, add weights wherever it goes. 
            % We will limit the weight per step.
            for ind = 1:length(active_agents)
                w = active_agents(ind);
                [inds, to_add] = assign_weights_basis(agent_directions(w),[istart(ind),jstart(ind)],[iend(ind),jend(ind)],phi_0,phi_grid_X,phi_grid_Y,weights,basis);
                % Update the weights.
                for k = 1 : size(inds,1)
                    if within_bounds_phi_mask(inds(k))
                        weights(inds(k)) = weights(inds(k)) + params.kappa*to_add(k);
                    end
                end 
            end
            
            % Update the recorded indices of the rounded points.
            agent_indices(active_agents,:) = endpoints;

            % Check for merge events via pairwise comparison.
            for w1 = active_agents(new_collisions)
                for w2 = active_agents(new_collisions' & active_agents>w1)
                    if all(agent_indices(w1,:) == agent_indices(w2,:),'all')
                        % A collision between agents has occured.
                        % Remove w1, the agent that didn't manage to draw.
                        is_active(w1) = false;
                        turned_off = w1;
                        % Record merge event.
                        merges = [merges; step_ind, w1, w2, turned_off];
                        % Remove the possibility of a collision with a trail for both agents.
                        % Note there is an edge case of a merge + collision happening in one go,
                        % but we ignore this here as it is very unlikely.
                        new_collisions([find(active_agents == w1);find(active_agents == w2)]) = false;
                    end
                end
            end

            % If agents are still active, if they were listed as a possible
            % collision, then they have collided with a trail. Deactive these
            % agents.
            collisions = [collisions; repmat(step_ind,sum(new_collisions),1), active_agents(new_collisions)', grid_X(endpoints(new_collisions,1))', grid_Y(endpoints(new_collisions,2))'];
            is_active(active_agents(new_collisions)) = false;

            % Allow active agents to split.
            active_agents = find(is_active);
            % If we might not have enough random choices, generate more.
            if rand_choice_counter + length(active_agents) > numRand
                [rand_choices, rand_direction, rand_choice_counter, rand_direction_counter] = generateRandoms(numRand,params);
            end
            % If the splitting counter has hit 0, split the agent.
            for w = active_agents
                if ~splitting_counters(w)
                    % Duplicate the coordinates.
                    N_agents = N_agents + 1;
                    agent_indices(N_agents,:) = agent_indices(w,:);
                    agent_positions(N_agents,:) = agent_positions(w,:);

                    % Record split event.
                    splits = [splits; step_ind, w, N_agents, agent_positions(N_agents,:)];

                    % Random +-.
                    sgn = rand_choices(rand_choice_counter); 
                    rand_choice_counter = rand_choice_counter + 1;

                    % Option 1.
                    % Generate the new directions, randomly picking which agent is turned which way.
                    agent_directions(N_agents) = agent_directions(w) + sgn * params.splitting_angle;
                    agent_directions(w) = agent_directions(w) - sgn * params.splitting_angle;

                    % Option 2.
                    % Let the original agent continue, and splice off the other.
                    % agent_directions(N_agents) = agent_directions(w) + sgn * params.splitting_angle;

                    % Move the newly spawned agent to a neighbouring free space in
                    % the collision grid in the direction of its motion. If none
                    % exists, leave the agent where it is, which may result in
                    % collision at the next timestep.
                    % Get the coordinates of the collision grid points around the current location.
                    nearby_indices = agent_indices(N_agents,:) + [-1,1;0,1;1,1;-1,0;1,0;-1,-1;0,-1;1,-1];
                    [i,j,lin_inds] = sanitise_indices(nearby_indices(:,1), nearby_indices(:,2), params.grid_X_num, params.grid_Y_num, within_bounds_mask);
                    % Check if positions are empty.
                    mask = ~trails(lin_inds);
                    if any(mask)
                        i = i(mask); j = j(mask);
                        % Choose a new space in the direction of the heading of the new
                        % agent, whilst also away from the heading of the parent.
                        directions = [grid_X(i)',grid_Y(j)'] - agent_positions(N_agents,:);
                        directions = directions ./ (sum(directions.^2,2).^0.5);
                        dot_prods = directions(:,1)*cos(agent_directions(N_agents)) + directions(:,2)*sin(agent_directions(N_agents));
                        mask = dot_prods > 0; i = i(mask); j = j(mask); directions = directions(mask,:);
                        if any(mask)
                            dot_prods = directions(:,1)*cos(agent_directions(w)) + directions(:,2)*sin(agent_directions(w));
                            [~,best_point] = min(abs(dot_prods));
                            agent_indices(N_agents,:) = [i(best_point), j(best_point)];
                            agent_positions(N_agents,:) = [grid_X(i(best_point)), grid_Y(j(best_point))];
                            % Fill in the location in trails.
                            trails = draw_agent(N_agents,[i(best_point), j(best_point)],trails);
                        end
                    end

                    % Generate new counters.
                    splitting_counters(N_agents) = 0;
                    % Activate the new agent.
                    is_active(N_agents) = true;
                end
            end

            % Decrement counters.
            splitting_counters = splitting_counters - 1;
            splitting_counters(splitting_counters<0) = params.splitting_period;

            active_agents = find(is_active);
            % Deactivate agents that reach the edges.
            for w = active_agents
                if ~within_bounds_mask(agent_indices(w,1),agent_indices(w,2))
                    is_active(w) = false;
                end
            end

            if ~any(active_agents)
                break
            end
        end

        weights_each_rep{rep_num} = weights;
        trails_each_rep{rep_num} = trails;
        event_log{rep_num} = struct('merges',merges,'splits',splits,'collisions',collisions);

        if params.num_repeats > 1
            textprogressbar(100*rep_num/params.num_repeats)
        end
    end

    % Package output.
    grids = {};
    grids.grid_X = grid_X;
    grids.grid_Y = grid_Y;
    grids.mesh_X = mesh_X;
    grids.mesh_Y = mesh_Y;
    grids.phi_mesh_X = phi_mesh_X;
    grids.phi_mesh_Y = phi_mesh_Y;
    grids.phi_grid_X = phi_grid_X;
    grids.phi_grid_Y = phi_grid_Y;
    grids.within_bounds_mask = within_bounds_mask;
    grids.within_bounds_phi_mask = within_bounds_phi_mask;
    grids.numbered_phi_grid = numbered_phi_grid;

    phi = {};
    phi.phi_0 = phi_0;
    phi.weights = weights_each_rep;

    trails = trails_each_rep;
end