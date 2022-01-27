function biases = compute_biases(phi, basis, grids, params)

    LR_biases = zeros(1,params.num_repeats+1);
    radial_biases = zeros(1,params.num_repeats+1);

    LR_biases(1) = measure_LR_bias(phi.phi_0);
    radial_biases(1) = measure_radial_bias(phi.phi_0, grids.phi_mesh_X, grids.phi_mesh_Y);
    for rep_ind = 1 : params.num_repeats
        lin_inds = 1 : numel(phi.phi_0);
        phi_current = reshape(eval_phi_basis(lin_inds, phi.phi_0(lin_inds(:)), basis, phi.weights{rep_ind}),size(phi.phi_0));

        LR_biases(rep_ind+1) = measure_LR_bias(phi_current);
        radial_biases(rep_ind+1) = measure_radial_bias(phi_current, grids.phi_mesh_X, grids.phi_mesh_Y);
    end

    biases = {};
    biases.LR_biases = LR_biases;
    biases.radial_biases = radial_biases;
end