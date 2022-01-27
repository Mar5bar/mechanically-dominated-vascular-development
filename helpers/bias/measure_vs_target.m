function d = metric_vs_target(candidate, target)
%% Compute a measure of difference between a candidate field phi and a target phi.
    d = nanmean(abs(cos(candidate - target)),'all');

end