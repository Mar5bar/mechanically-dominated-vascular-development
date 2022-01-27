function counters = init_splitting_counters(splitting_period,N_walkers)
% INIT_SPLITTING_COUNTERS assigns a random initial counter value from 0,...,period to each walker.
    counters = randi(splitting_period+1,N_walkers,1) - 1;
end