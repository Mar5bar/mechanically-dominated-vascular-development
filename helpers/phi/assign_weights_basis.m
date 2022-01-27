function [lin_inds, to_add] = assign_weights_basis(direction,I1,I2,phi_0,phi_X_grid,phi_Y_grid,weights,basis)
	d = I2 - I1;
    num_steps = sum(abs(d));
    if num_steps == 0
        inds = I1;
    else
        inds = zeros(num_steps,2);
        d = d / num_steps;
        inds(1,:) = I1 + d;
        for i = 2 : num_steps
            inds(i,:) = inds(i-1,:) + d;
        end
        inds = round(inds);
    end
    % Compute the current phi at these points.
    lin_inds = inds(:,1) + (inds(:,2)-1) * size(phi_X_grid,2);
    phi = eval_phi_basis(lin_inds,phi_0(lin_inds),basis,weights);
    % Compute the weight sign and strength.
    v = [cos(direction),sin(direction)];
    w = [cos(phi),sin(phi)];
    dot_prod = w*v';

    % This mask is 1 if vectors are pointing in the same direction, and -1
    % if they are pointing away from each other and one needs reversing.
    mask = dot_prod>0; mask = 2*mask-1;
    % If facing the wrong way, negate v.
    swapped_vs = mask.*v;
    sgn = -sign(diff(unwrap([atan2(swapped_vs(:,2),swapped_vs(:,1)),phi],[],2),[],2));
    mag = real(acos(abs(dot_prod)));

    % The magnitude of the weight should be such that it corresponds to
    % enacting the change in phi so as to perfectly align phi with the
    % direction. This can then be modulated by multiplication by 0<=kappa<=1.
    to_add = mag .* sgn;
end