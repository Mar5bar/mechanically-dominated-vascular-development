function [trails, possible_collisions] = draw_agents(labels,indices,trails)
% Do not check any bounds.
    if numel(labels) == 1
        labels = labels * ones(size(indices,1));
    end
    if nargout == 2
        possible_collisions = logical(zeros(size(indices,1),1));
        for i = 1 : size(indices,1)
            [trails, possible_collisions(i)] = draw_agent(labels(i),indices(i,:),trails);
        end
    else
        for i = 1 : size(indices,1)
            trails = draw_agent(labels(i),indices(i,:),trails);
        end
    end
end