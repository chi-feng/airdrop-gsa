function threshold = survival(vf, mu, sigma)
% computes the survival probability for an impact velocity of vf given 
% parameters mu and sigma. Note the 50% surivival rate is at vf = exp(mu)

threshold = 0.5 - 0.5 * erf((log(vf) - mu) / (sqrt(2) * sigma));

end