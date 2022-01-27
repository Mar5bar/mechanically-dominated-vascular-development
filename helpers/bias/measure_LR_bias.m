function bias = measure_LR_bias(phi)
    % Compute a measure of LR bias. For uniformly distributed phi on
    % [0,2pi], this value should be 2/pi.
    bias = nanmean(abs(cos(phi(:))));
end