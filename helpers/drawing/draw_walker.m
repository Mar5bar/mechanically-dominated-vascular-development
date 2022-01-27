function [trails, overwriting] = draw_walker(label,I,trails,allowed)
% DRAW_WALKER draws a walker on trails at (i,j).
% We're not going to check any bounds, but we will check for overwrites/collisions.
    ind = I(1) + (I(2)-1)*size(trails,1);
    if nargout == 2
        % Flag up if we're about to overwrite anything made by another walker.
        test_val = trails(ind);
        if nargin < 4
            overwriting = test_val ~=0;
        else
            % Allow collisions with 'allowed', and overwrite.
            overwriting = test_val ~=0 & test_val~=allowed;
        end
        % Only draw if we wouldn't be overwriting anything.
        if ~overwriting
            trails(ind) = label;
        end
    else
        trails(ind) = label;
    end
end
