%% Parameters.
params = {};
params.N_steps = 1000;  % Number of timesteps.
params.N_walkers = 20;  % Number of initial walkers.
params.speed = 0.01;    % Linear speed of the walkers, here fixed. Should be such that walkers do not attempt to draw over themselves.
params.k = 0.8;         % The strength of alignment of walkers to the phi field, in [0,1].
params.kappa = 0.8;       % Effectiveness of remodelling, between 0 and 1. 1 corresponds to perfect remodelling.
params.lambda = 10;   % Spatial parameter of basis function.
params.noise_variance = 1/100;  % variance of normal distribution for rand_direction. Must be scaled with timestep.
params.splitting_period = 25;   % Number of steps needed until a split occurs.
params.splitting_angle = pi/6;  % Half-angle of split walkers, in radians,
params.init_phi = 'horizontal';     % Initial phi, one of 'radial', 'random', 'horizontal', 'vertical', or 'spiral'.
params.num_repeats = 1;         % Number of times to repeat simulation.
params.reset_phi = true;        % Initialise weights to zero. Set to false to examine regrowth.
params.degrade_phi = false;     % Degrade weights by 1/2 at the beginning of the simulation, modelling the decay of memory.
% Grids
params.grid_X_lims = [-1,1];
params.grid_X_num = 1000;
params.grid_Y_lims = [-1,1]; 
params.grid_Y_num = 1000; 
params.phi_grid_X_num = 100;
params.phi_grid_Y_num = 100;
% Domain
params.innerRadius = 0.1;
params.outerRadius = 1;

if params.num_repeats > 1 && params.reset_phi
    disp("NOTE: Phi is being reset every repeat. Set reset_phi = false to preserve phi.")
end

addpath(genpath('./'))

tic 
[trails, phi, grids, event_log, basis] = runSim(params);
toc

% Evaluate biases
biases = compute_biases(phi, basis, grids, params);

%% Save as most recent output.
save('output/recent_output.mat','trails','phi','grids','event_log','biases','params')

%% Plot network after final repeat.
plot_network(trails, phi, grids, params, params.num_repeats)
drawnow