function indices = gen_circular_mask_indices(radius)
% GEN_CIRCULAR_MASK_INDICES creates the index offsets necessary for drawing
% a circle with specified radius about a point.
    xs = -radius:radius; ys = xs;
    [xgrid, ygrid] = meshgrid(xs,ys);
    mask = xgrid.^2 + ygrid.^2 <= radius^2;
    xs = xgrid(mask); ys = ygrid(mask);
    indices = [xs(:), ys(:)];
end