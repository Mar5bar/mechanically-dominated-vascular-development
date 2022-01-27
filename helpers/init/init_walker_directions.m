function thetas = init_walker_directions(positions)
% INIT_WALKER_DIRECTIONS generates and returns the raw directions of the
% initial walkers.
	thetas = atan2(positions(:,2),positions(:,1));
end
