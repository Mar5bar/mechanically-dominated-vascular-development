function bias = measure_radial_bias(phi,X,Y)
    % Compute a measure of LR bias. For radial phi, this value should be 1.
    % For random phi, this should be 2/pi.
    phi_radial = atan2(Y,X);
    bias = nanmean(abs(cos(phi(:) - phi_radial(:))));
end