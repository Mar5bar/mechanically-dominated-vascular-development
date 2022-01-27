function [i,j,lin_inds] = sanitise_indices(i,j,X_num,Y_num,mask)
% SANITISE_INDICES checks that the indices i,j are within the grid of
% dimensions (X_num,Y_num) and are such that mask(i,j) = true.
    lin_inds = i + (j-1) * X_num;
    valid_mask = ~(i <= 0 | j <= 0 | i > X_num | j > Y_num);
    if any(valid_mask)
        valid_mask(valid_mask) = mask(lin_inds(valid_mask));
    end
    i = i(valid_mask);
    j = j(valid_mask);
    lin_inds = lin_inds(valid_mask);
end