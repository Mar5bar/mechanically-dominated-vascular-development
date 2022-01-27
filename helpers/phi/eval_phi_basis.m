function phi = eval_phi_basis(lin_inds,phi_0,basis,weights)
    lin_inds = lin_inds(:);
    phi_0 = phi_0(:);
    phi = (weights*basis(:,lin_inds))' + phi_0;
end