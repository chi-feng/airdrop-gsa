function dudt = payload_dudt(t, u, m, r, Cd, wx, tfree, topen)

rp      = 0.5;  % radius of parachute in meters
Cdp     = 1.75; % parachute coefficient of drag

rho  = 1.22;    % atmospheric density
g    = 9.8;     % gravity
mu   = 1.81e-5; % viscosity of air 
Fmax = 300;     % maximum force on parachute (N)

% unwrap u
x    = u(1);
y    = u(2);
vx   = u(3);
vy   = u(4);
detached = u(5); % 0 if intact, > 0 if torn

% calculate relative velocity 
v = sqrt((vx - wx)^2 + vy^2);

% determine parachute opening status
if t < tfree || detached > 0
    parachute = 0;
elseif t > (tfree + topen)
    parachute = 1;
else
    parachute = min(1, (t - tfree) / topen);
end

% compute Reynolds number
Re = rho * v * r / mu;

% compute atmospheric drag due to payload and parachute
if abs(Re) > 0
    Cd_eff = 24 / Re + 6 / (1 + sqrt(Re)) + Cd; 
    Fd = -0.5 * rho * pi * r^2 * v^2 * Cd_eff;
    Cdp_eff = 24 / Re + 6 / (1 + sqrt(Re)) + Cdp;
    Fdp = -0.5 * rho * pi * rp^2 * v^2 * Cdp_eff * parachute; 
else
    Fd = 0;
    Fdp = 0;
end

% approximate effective drag as sum
ax = -abs(Fd + Fdp) / m * (vx - wx) / v;
ay = -g - abs(Fd + Fdp) / m * vy / v;

dudt = [ vx, vy, ax, ay, abs(Fdp) > Fmax]';

end