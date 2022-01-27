function [trails, collisions, endpoints] = draw_paths(labels,starts,ends,trails)
    collisions = logical(zeros(size(starts,1),1));
    endpoints = zeros(size(starts,1),2);
    for i = 1 : size(starts,1)
        [trails, collisions(i), endpoints(i,:)] = draw_path(labels(i),starts(i,:),ends(i,:),trails);
    end
end