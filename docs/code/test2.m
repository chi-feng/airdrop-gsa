clear all
close all
v = linspace(0.01,20,100);
mu = 2.3;
sigma = 0.15;
thresh = 0.5 - 0.5 * erf((log(v) - mu) / (sqrt(2) * sigma));
plot(v,thresh, 'LineWidth',2)
ylabel('Probability of surviving impact')
xlabel('Impact velocity (m/s)')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 4.5 3];
print('../figures/survival','-depsc')