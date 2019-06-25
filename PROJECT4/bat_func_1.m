function val = bat_func_1(xv)

syms x y

eq1 = ((x/7)^2*sqrt(abs(abs(x)-3)/(abs(x)-3))+(y/3)^2*sqrt(abs(y+3/7*sqrt(33))/(y+3/7*sqrt(33)))-1);
%eq2 = (abs(x/2)-((3*sqrt(33)-7)/112)*x^2-3+sqrt(1-(abs(abs(x)-2)-1)^2));
eq3 = (9*sqrt(abs((abs(x)-1)*(abs(x)-.75))/((1-abs(x))*(abs(x)-.75)))-8*abs(x));
eq4 = (3*abs(x)+.75*sqrt(abs((abs(x)-.75)*(abs(x)-.5))/((.75-abs(x))*(abs(x)-.5))));
eq5 = (2.25*sqrt(abs((x-.5)*(x+.5))/((.5-x)*(.5+x))));
eq6 = (6*sqrt(10)/7+(1.5-.5*abs(x))*sqrt(abs(abs(x)-1)/(abs(x)-1))-(6*sqrt(10)/14)*sqrt(4-(abs(x)-1)^2));


% eq2 is tricky to evaluate, so it's quickest to represent it using a
% spline. 
xf = [3  ,  4,  5.2, 6.7, 6.9, 7];
yf = [2.7,2.5, 2.05, .95, 0.5, 0];

for n = 1:length(xv)
    
    xd = xv(n);
    
    xd = abs(xd);
    
    % there's a singularity at 1
    if xd == 1
        xd = 1.00000000001;
    end
    
    % and one at 0.5
    if xd == 0.5
        xd = 0.49999999999;
    end
    
    % and one at 0.75
    if xd == 0.75
        xd = 0.749999999999;
        
    end
    
    if      xd <= 0.5
        
        val(n) = subs(eq5,xd);
        
    elseif  xd > 0.5 && xd <= 0.75
        
        val(n) = subs(eq4,xd);
        
    elseif  xd > 0.75 && xd <= 1.0

        val(n) = subs(eq3,xd);
        
    elseif  xd > 1.0 && xd <= 3.0
        
        val(n) = subs(eq6,xd);
        
    elseif  xd > 3.0
        
        % spline value at xd
        val(n) = pchip(xf,yf,xd);
        
    end
    
end


