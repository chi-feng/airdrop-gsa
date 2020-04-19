function [position, isterminal, direction] = ground_intersection(t, u)
% checks if payload hit the ground

position   = u(2); % y(t)
isterminal = 1;    % halt integration
direction  = 0;    % approach from either direction

end