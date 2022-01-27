function [trails, flag, endpoint] = draw_path_fast(label,I1,I2,trails)
% DRAW_PATH draws a path in trails connecting I1 to I2.
    d = I2 - I1;
    num_steps = sum(abs(d));
    if num_steps == 0
        flag = false;
        return
    end
    inds = zeros(num_steps,2);
    d = d / num_steps;
    inds(1,:) = I1 + d;
    for i = 2 : num_steps
        inds(i,:) = inds(i-1,:) + d;
    end
    inds = round(inds);
    % Remove duplicates.
    inds([false;all(~diff(inds),2)],:) = [];
    % Remove the start point.
    inds(all(inds==I1,2),:) = [];
    lin_inds = inds(:,1) + (inds(:,2)-1)*size(trails,1);
    can_draw = trails(lin_inds) == 0;
    flag = any(~can_draw);
    if flag
        last_ind = find(~can_draw,1,'first')-1;
        if last_ind == 0
            endpoint = I1;
        else
            endpoint = inds(last_ind,:);
        end
        can_draw = can_draw(1:last_ind);
    else
        endpoint = I2;
    end
    trails(lin_inds(can_draw)) = label;
end
