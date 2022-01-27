function counters = init_splitting_counters(splitting_period,N_agents)
% INIT_SPLITTING_COUNTERS assigns a random initial counter value from 0,...,period to each agent.
    counters = randi(splitting_period+1,N_agents,1) - 1;
end