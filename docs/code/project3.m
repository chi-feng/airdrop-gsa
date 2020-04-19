% nominal values of parameters
m       = 6;    % mass of payload and parachute assembly kilograms
r       = 0.1;  % radius of payload in meters
Cd      = 0.5;  % payload coefficient of drag
wx      = 0;    % horizontal wind speed in m/s
tfree   = 9;    % time before parachute opens
topen   = 5;    % time it takes for parachute to open completely

% initial conditions x(0), y(0), vx(0), vy(0)
u0 = [-340, 500, 50, 0];

% compute a single trajectory
[t, u] = payload_sim(u0, m, r, Cd, wx, tfree, topen);

% plot the trajectory
close all;
figure(1);
plot(u(:,1), u(:,2), '-')
xlabel('horizontal displacement (m)');
ylabel('altitude (m)');

% plot the velocity over time
figure(2);
plot(t, sqrt(u(:,3).^2 + u(:,4).^2),'-')
xlabel('time after release (s)');
ylabel('velocity (m/s)');

%% Monte Carlo Simulation

% Monte Carlo loop
N = 1000;

% initialize input samples

x  = zeros(N,1);
y  = zeros(N,1);
vx = zeros(N,1);
vy = zeros(N,1);
u0 = zeros(N,4);

% initialize output samples
xf     = zeros(N,1); % landing sites (m)
vf     = zeros(N,1); % impact velocities (m/s)
intact = zeros(N,1); % payload intact after impact (0, 1)

for i=1:N % consider using parfor, especially for large N
    
    % Draw samples from the input distributions
    
    % initial conditions
    x(i)     = ; 
    y(i)     = ;
    vx(i)    = ;
    vy(i)    = ;
    u0(i,:)  = [x(i), y(i), vx(i), vy(i)];
    
    % random inputs
    m(i)     = ;
    r(i)     = ;
    Cd(i)    = ;
    tfree(i) = ;
    topen(i) = ; 
    wx(i)    = ;
    
    % simulate trajectory
    [t, u] = payload_sim(u0(i,:), m(i), r(i), Cd(i), wx(i), tfree(i), topen(i));
    
    % get the impact location and impact velocity
    xf(i) = u(end,1);
    vf(i) = sqrt(u(end,3)^2 + u(end,4)^2);
    
    % check if payload survived impact
    if rand() < survival(vf(i), mu, sigma)
        intact(i) = 1;
    end
    
end

%% Monte Carlo Analysis

% Below is starting point for conducting your Monte Carlo analysis

inside = find(abs(xf) < 50);
outside = find(abs(xf) > 50);

survived = find(intact == 1);
destroyed = find(intact == 0);

survived_inside = intersect(survived, inside);

% probability of success (outcome A)
p_success = length(survived_inside) / N;

% plot histogram of landing sites of intact payloads
histogram(xf(survived),'BinWidth',5);
