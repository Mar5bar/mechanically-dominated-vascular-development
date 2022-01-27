function phi = init_phi(X,Y,mask,init_phi_str)
% INIT_THETA will generate the initial vector phi_0 containing orientation information.
    phi = zeros(size(X));

    switch init_phi_str
    case 'radial'
        % Radial
        ang = atan2(Y,X);
        phi(mask) = ang(mask);
    case 'random'
        % Random.
        phi(mask) = rand(nnz(mask),1)*pi;
    case 'horizontal'
        % Horizontal.
        phi(mask) = 0;
    case 'vertical'
        % Vertical.
        phi(mask) = pi/2;
    case 'spiral'
        % Spiral.
        phi = atan2(-X-Y,Y-X);
    end
    phi(~mask) = NaN;
end