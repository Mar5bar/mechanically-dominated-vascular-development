function [rand_choices, rand_direction, rand_choice_counter, rand_direction_counter] = generateRandoms(numRand,params)
%% Generate the required arrays of random numbers and choices.
        rand_choices = 2*round(rand(1,numRand))-1;
        rand_choice_counter = 1;
        rand_direction = normrnd(0,sqrt(params.noise_variance),numRand,1); % Precompute choices.
        rand_direction_counter = 1;
end