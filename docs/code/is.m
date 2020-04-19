openclear all;
close all;

nsamp = 1000;
replicates = 1000;

f = @(x) exp(-(x-pi/2).^2/0.01);

% f = @(x) (abs(x-pi/2)<0.1)*1;

p = @(x) tripdf(x,-pi,0,pi);

mu = pi/2;
sigma = 0.1;
q = @(x) 1/(sqrt(2*pi*sigma^2))*exp(-(x-mu).^2/(2*sigma^2));

fill_between_lines = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C );

x_coarse = linspace(-pi,pi,101);
x_fine = linspace(pi/4,3*pi/4,400);

x = union(x_coarse,x_fine,'sorted');

exact = trapz(x,f(x).*p(x));
%exact = 0.01994711402007164;

figure;
plot(x,f(x),'Color',[0.1 0.3 0.8]) 
hold on
plot(x,p(x),'Color',[0.8 0.3 0.1])
plot(x,f(x).*p(x),'Color',[0.6 0.1 0.7])
plot(x,q(x),'Color',[0.6 0.6 0.6])
legend('$f(x)$','$p(x)$','$f(x)p(x)$','$q(x)$','Location','NorthWest')

h = fill_between_lines(x,f(x).*p(x),x*0,[0.6 0.1 0.7]);
set(h,'edgealpha',0)
set(h,'facealpha',0.5)
hold off

xticks([-pi,0,pi/2,pi])
xticklabels({'-\pi','0','\pi/2','\pi'})
xlim([-pi,pi])
ylim([0,1])

estimates = zeros(replicates,nsamp);
summands = [ ];
for i = 1:replicates
    x = trirnd(-pi,0,pi,nsamp);
    if i <= 100
        summands = [summands; f(x)];
    end
    estimates(i,:) = cumsum(f(x))'./(1:nsamp);
end

figure;
p95 = fill_between_lines(1:nsamp,prctile(estimates,95,1),prctile(estimates,5,1),[0.8 0.8 0.8]);
set(p95,'edgealpha',0)
hold on
plot(1:nsamp,estimates(1:5,:)')
h = plot(1:nsamp,mean(estimates,1),'LineWidth',1.5,'Color',[0 0 0]);
hold off
ylim([0,0.05])
legend([h,p95],{'Estimator mean','5th and 95th percentile'},'Location','NorthEast');

figure;
histogram(summands,'Normalization','probability')
xlabel('Importance sampling estimator summand value')
ylabel('Probability');

nsamps = floor(logspace(2,4,5));
mse = zeros(length(nsamps),1);
estimates = zeros(length(nsamps),replicates);
for j = 1:length(nsamps)
    for i = 1:replicates
        x = trirnd(-pi,0,pi,nsamps(j));
        estimates(j,i) = sum(f(x))/nsamps(j);
    end 
    mse(j) = sum((estimates(j,:) - exact).^2)/replicates;
end

figure;
loglog(nsamps,(mse),'.-');
ylabel('MSE');
xlabel('N');
hold on

nsamps = floor(logspace(2,4,5));
mse = zeros(length(nsamps),1);
estimates = zeros(length(nsamps), replicates);
for j = 1:length(nsamps)
    for i = 1:replicates
        x = randn(nsamps(j),1)*sigma + mu;
        w = p(x)'./q(x);
        estimates(j,i) = sum(w.*f(x))/nsamps(j);
    end 
    mse(j) = sum((estimates(j,:) - exact).^2)/replicates;
end
loglog(nsamps,(mse),'.-');
hold off
legend('Direct Monte Carlo','Importance Sampling');


estimates = zeros(replicates,nsamp);
summands = [ ];
for i = 1:replicates
    x = randn(nsamp,1)*sigma + mu;
    w = p(x)'./q(x);
    if i <= 100
        summands = [summands; w.*f(x)];
    end
    estimates(i,:) = cumsum(w.*f(x))'./(1:nsamp);
end

figure;
p95 = fill_between_lines(1:nsamp,prctile(estimates,95,1),prctile(estimates,5,1),[0.8 0.8 0.8]);
set(p95,'edgealpha',0)
hold on
plot(1:nsamp,estimates(1:5,:)')
h = plot(1:nsamp,mean(estimates,1),'LineWidth',1.5,'Color',[0 0 0]);
hold off
ylim([0,0.05])
legend([h,p95],{'Estimator mean','5th and 95th percentile'},'Location','NorthEast');

figure;
histogram(summands,'Normalization','probability')
xlabel('Importance sampling estimator summand value')
ylabel('Probability');



function p = tripdf(x,a,c,b)
    p = zeros(1,length(x));
    left = find((a<=x)&(x<c));
    right = find((c<=x)&(x<=b));
    p(left) = 2*(x(left)-a)/(b-a)/(c-a);
    p(right) = 2*(b-x(right))/(b-a)/(b-c);
end

function x = trirnd(a,c,b,n)
    if nargin == 3
        n = 1;
    end
    x = zeros(n,1);
    u = rand(n,1);
    F = (c-a)/(b-a);
    smaller = find(u<F);
    larger = find(u>=F);
    x(smaller) = a+sqrt(u(smaller)*(b-a)*(c-a));
    x(larger) = b-sqrt((1-u(larger))*(b-a)*(b-c));
end

