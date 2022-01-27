function new_theta = gen_agent_direction(agent,theta,phi,noise,k)
% GEN_WALKER_DIRECTION generates the agent direction new_theta given the old direction, theta,
% and the orientation field phi.
    if isnan(phi)
        d_angle = 0;
    else
        v = [cos(theta),sin(theta)];
        w = [cos(phi),sin(phi)];
        dot_prod = w*v';

        if dot_prod < 0
            v = -v;
        end

        sgn = sign(diff(unwrap([atan2(v(:,2),v(:,1)),phi],[],2),[],2));
        mag = real(acos(abs(dot_prod)));
        d_angle = sgn * mag;
    end
    
    new_theta = theta + k*d_angle + noise;
end