x = 1:100;
y = sin(x);
e = std(y)*ones(size(x));
errorbar(x,y,e)
set(gca,'yscale','log')