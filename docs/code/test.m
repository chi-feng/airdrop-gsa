%clear all;
make_plot = 0;

mode = 2;

if mode == 2
    load biasing_dist;
end

if make_plot

% nominal values of parameters
m       = 5;    % mass of payload and parachute assembly kilograms
r       = 0.1;  % radius of payload in meters
Cd      = 0.5;  % payload coefficient of drag
wx      = 0;    % horizontal wind speed in m/s
tfree   = 9;      % time before parachute opens
topen   = 5;    % time it takes for parachute to open completely

% initial conditions x(0), y(0), vx(0), vy(0)
u0 = [-326.55, 500, 50, 0];

% compute a single trajectory
[t, u] = payload_sim(u0, m, r, Cd, wx, tfree, topen);

close all;

figure(1);
plot(u(:,1), u(:,2), '-')
xlabel('horizontal displacement (m)');
ylabel('altitude (m)');
[~, idx_topen] = min(abs(t-(tfree+topen)));
[~, idx_tfree] = min(abs(t-(tfree)));
text(u(idx_tfree,1), u(idx_tfree,2), 't_{free}\rightarrow', 'HorizontalAlignment', 'right')
text(u(idx_topen,1), u(idx_topen,2), 't_{free} + t_{open}\rightarrow', 'HorizontalAlignment', 'right')

ylim([0,500])
xlim([-400,100])
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 4.5 3];
% print('../figures/trajectory','-depsc')

figure(2);
ax = axes;
plot(t, sqrt(u(:,3).^2 + u(:,4).^2),'-')
hold on
xlim([0,26]);
ylim([0,70]);
line([tfree tfree],[0,70],'Color',[0.5 0.5 0.5])
line([tfree+topen tfree+topen],[0,70],'Color',[0.5 0.5 0.5])
text(tfree, 5, 't_{free}\rightarrow', 'HorizontalAlignment', 'right')
text(tfree + topen, 4, '\leftarrow t_{free} + t_{open}', 'HorizontalAlignment', 'left')
xlabel('time after release (s)');
ylabel('velocity (m/s)');
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 4.5 3];
% print('../figures/velocity','-depsc')

end

N = 10000;

u0 = zeros(N,4);
xf = zeros(N,1);
vf = zeros(N,1);
pf = zeros(N,1);
intact = ones(N,1);

x       = zeros(N,1);
y       = zeros(N,1);
vx      = zeros(N,1);
vy      = zeros(N,1);
m       = zeros(N,1);
r       = zeros(N,1);
Cd      = zeros(N,1);
tfree   = zeros(N,1);
topen   = zeros(N,1);
wx      = zeros(N,1);

if mode == 2
    load('biasing_dist');
    inputs = mvnrnd(biasing_dist.mean,biasing_dist.cov,N);
    q = mvnpdf(inputs,biasing_dist.mean,biasing_dist.cov);
end

parfor i = 1:N
    
    if mode == 1 % standard
        x(i)       = randn() * 3 - 340;
        y(i)       = randn() * 3 + 500;
        vx(i)      = randn() * 0.5 + 50;
        vy(i)      = randn() * 0.5;
        m(i)       = trirand(4,6,9);
        r(i)       = trirand(0.09,0.1,0.11);
        Cd(i)      = trirand(0.45,0.5,0.6);
        tfree(i)   = trirand(8.95,9.0,9.05);
        topen(i)   = logrand(1.65, 0.15); 
        wx(i)      = 2 * randn();
    end
    
    if mode == 2 % importance sampling from file
        x(i)       = inputs(i,7);
        y(i)       = inputs(i,8);
        vx(i)      = inputs(i,9);
        vy(i)      = inputs(i,10);
        m(i)       = inputs(i,1);
        r(i)       = inputs(i,2);
        Cd(i)      = inputs(i,3);
        tfree(i)   = inputs(i,5);
        topen(i)   = inputs(i,6); 
        wx(i)      = inputs(i,4);
    end
    
    if mode == 6 % importance sampling
        m(i)       = logrand(1.6, 0.001); 
        r(i)       = trirand(0.09,0.1,0.11);
        Cd(i)      = trirand(0.45,0.5,0.6);
        tfree(i)   = trirand(8.99,9.0,9.01);
        topen(i)   = logrand(1.65, 0.01); 
        wx(i)      = 0.01 * randn();
        x(i)       = randn() * 0.01 - 325;
        y(i)       = randn() * 0.01 + 500;
        vx(i)      = randn() * 0.01 + 50;
        vy(i)      = randn() * 0.01;
    end
    
    u0(i,:)    = [x(i), y(i), vx(i), vy(i)];
    
    [t, u] = payload_sim(u0(i,:), m(i), r(i), Cd(i), wx(i), tfree(i), topen(i));
    xf(i) = u(end,1);
    vf(i) = sqrt(u(end,3)^2+u(end,4)^2);
    pf(i) = u(end,5);
    
    mu = 2.3;
    sigma = 0.15;
    thresh = 0.5 + 0.5 * erf((log(vf(i)) - mu) / (sqrt(2) * sigma));
    if rand() < thresh
        intact(i) = 0;
    end
end


%%
close all 


inside = find(abs(xf) < 1);
outside = find(abs(xf) > 1);
survived = find(intact == 1);
destroyed = find(intact == 0);

survived_inside = intersect(inside, survived);
survived_outside = intersect(outside, survived);
destroyed_inside = intersect(inside, destroyed);
destroyed_outside = intersect(outside, destroyed);

success = survived_inside;
% success = survived_inside;

if mode == 2
    is_samples = zeros(N,1);
    for j = 1:length(success)
        i = success(j);
        nominal = [
            log(normpdf(x(i),-340,3))
            log(normpdf(y(i),500,3))
            log(normpdf(vx(i),50,0.5))
            log(normpdf(vy(i),0,0.5))
            log(tripdf(m(i),4,6,9))
            log(tripdf(r(i),0.09,0.1,0.11))
            log(tripdf(Cd(i),0.45,0.5,0.6))
            log(tripdf(tfree(i),8.95,9.0,9.05))
            log(lnormpdf(topen(i),1.65,0.15))
            log(normpdf(wx(i),0,2))
        ]';
        is_samples(i) = exp(sum(nominal) - log(q(i)));
    end
else
    mc_samples = zeros(N,1);
    mc_samples(success)=1;
end

figure
subplot(5,2,1); plot_input_histogram(m(:), m(success)); xlabel('m');
subplot(5,2,2); plot_input_histogram(r(:), r(success)); xlabel('r');
subplot(5,2,3); plot_input_histogram(Cd(:), Cd(success)); xlabel('Cd');
subplot(5,2,4); plot_input_histogram(wx(:), wx(success)); xlabel('wx');
subplot(5,2,5); plot_input_histogram(tfree(:), tfree(success)); xlabel('tfree');
subplot(5,2,6); plot_input_histogram(topen(:), topen(success)); xlabel('topen');
subplot(5,2,7); plot_input_histogram(u0(:,1), u0(success,1)); xlabel('x(0)'); 
subplot(5,2,8); plot_input_histogram(u0(:,2), u0(success,2)); xlabel('y(0)');
subplot(5,2,9); plot_input_histogram(u0(:,3), u0(success,3)); xlabel('vx(0)');
subplot(5,2,10); plot_input_histogram(u0(:,4), u0(success,4)); xlabel('vy(0)');

p_survived_inside = length(survived_inside) / N
p_survived_outside = length(survived_outside) / N
p_destroyed_inside = length(destroyed_inside) / N
p_destroyed_outside = length(destroyed_outside) / N

figure
hold on
histogram(vf(survived),'BinWidth',0.5)
histogram(vf(destroyed),'BinWidth',0.5)
xlim([0,20])
legend('survived','destroyed','Location','NorthWest')

figure
histogram(xf(survived),'BinWidth',1)
hold on
histogram(xf(destroyed),'BinWidth',1)
legend('survived','destroyed','Location','NorthWest')

%%

inputs = [m r Cd wx tfree topen u0];

rho = tril(corr(inputs(success,:)));
figure(4)
%surf(rho)
plotmatrix(inputs(success,:))
xticks(1:10)
yticks(1:10)
yticklabels({'m' ,'r' ,'Cd' ,'wx' ,'tfree' ,'topen' ,'x','y','vx','vy'})
xticklabels({'m' ,'r' ,'Cd' ,'wx' ,'tfree' ,'topen' ,'x','y','vx','vy'})

mu = mean(inputs(success,:));
sigma = cov(inputs(success,:));
save('q_mean','mu');
save('q_cov','sigma');

%%
figure(5)
plot((1:N)',cumsum(mc_samples) ./ (1:N)','LineWidth',1.5)
hold on
plot((1:N)',cumsum(is_samples) ./ (1:N)','LineWidth',1.5)
hold off
ylim([0,0.1])
legend('Monte Carlo', 'Importance Sampling')
%%
% 
% figure(4)
% % plot(wx,xf,'.')
% subplot(1,2,1)
% histogram2(wx,m,'FaceColor','flat','ShowEmptyBins','on');
% 
% subplot(1,2,2)
% histogram2(wx(success),m(success),'FaceColor','flat','ShowEmptyBins','on');
% %%
% % 
% inputs = [m' r' Cd' wx' tfree' topen' u0];
% rho = corr(inputs);
% figure(5)
% subplot(1,2,1)
% bar3(rho)
% xticks(1:10)
% yticks(1:10)
% yticklabels({'m' ,'r' ,'Cd' ,'wx' ,'tfree' ,'topen' ,'x','y','vx','vy'})
% xticklabels({'m' ,'r' ,'Cd' ,'wx' ,'tfree' ,'topen' ,'x','y','vx','vy'})
% subplot(1,2,2)
% rho_success = corr(inputs(success,:));
% bar3(rho_success)
% xticks(1:10)
% yticks(1:10)
% yticklabels({'m' ,'r' ,'Cd' ,'wx' ,'tfree' ,'topen' ,'x','y','vx','vy'})
% xticklabels({'m' ,'r' ,'Cd' ,'wx' ,'tfree' ,'topen' ,'x','y','vx','vy'})