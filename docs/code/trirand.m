function x = trirand(a, c, b)
% draw a sample from a triangular distribution with min=a, mode=b, max=c

fc = (c - a) / (b - a);
u = rand();
if u < fc
    x = a + sqrt(u * (b - a) * (c - a));
else
    x = b - sqrt((1-u)*(b-a)*(b-c));
end

end