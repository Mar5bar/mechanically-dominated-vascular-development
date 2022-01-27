function thetas = init_agent_directions(positions)
% INIT_WALKER_DIRECTIONS generates and returns the raw directions of the
% initial agents.
	thetas = atan2(positions(:,2),positions(:,1));
end
