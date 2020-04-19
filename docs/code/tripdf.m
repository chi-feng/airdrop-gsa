% compute pdf of triangular distribution
function p = tripdf(x,a,c,b)
    
    p = zeros(1,length(x));
    
    left = find((a<=x)&(x<c));
    right = find((c<=x)&(x<=b));
    
    p(left) = 2*(x(left)-a)/(b-a)/(c-a);
    p(right) = 2*(b-x(right))/(b-a)/(b-c);
    
end