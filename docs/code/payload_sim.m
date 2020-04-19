function [t, u] = payload_sim(u0, m, r, Cd, wx, tfree, topen)
% compute deterministic payload trajectory given input parameters

u0 = [u0, 0];
tspan = [0, 120];
options = odeset('Events', @ground_intersection, 'RelTol', 1e-3, 'AbsTol', 1e-3);
[t, u] = ode45(@(t, u) payload_dudt(t, u, m, r, Cd, wx, tfree, topen), tspan, u0, options);
u(:,5) = sign(u(:,5)); % threshold parachute_detached

end