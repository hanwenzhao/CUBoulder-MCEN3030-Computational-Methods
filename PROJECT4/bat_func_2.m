function val = bat_func_2(xv)

syms x y

eq2 = (abs(x/2)-((3*sqrt(33)-7)/112)*x^2-3+sqrt(1-(abs(abs(x)-2)-1)^2));

% eq2 is tricky to evaluate, so it's quickest to represent it using a
% spline. 
xf = [3  ,  4,  5.2, 6.7, 6.9, 7];
yf = [2.7,2.5, 2.05, .95, 0.5, 0];

for n = 1:length(xv)
    
    xd = xv(n);
    
    xd = abs(xd);
    
    
    if      xd <= 4
        
        val(n) = subs(eq2,xd);
        
    else
        
        %spline value at xd
        val(n) = -pchip(xf,yf,xd);
        
    end
    
end


