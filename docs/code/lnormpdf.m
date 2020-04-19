function f=lnormpdf(x,mu,sigma)
    f = 1/(x*sigma*sqrt(2*pi))*exp(-(log(x)-mu)^2/(2*sigma^2));
end