% MCEN 3030
% MEID: 650-703
% PROJECT 4
% HANWEN ZHAO
function [] = cp4_650703() % main function
clc
%% Problem 1
a = -7; % assign the boundary value
b = 7;
h = 0.1; % initial h

% calculate the first function
I_old1 = 0;
times1 = 0;
h1 = h;
fprintf('Evaluating bat_function_1.\n')
while 1
    times1 = times1 + 1; % use times to record the number of times to evaluate
    fprintf('Attempt %d.\n', times1)
    I_new1 = run_simp(@bat_func_1, a, b, h1); % call simpson's 1/3 method
    fprintf('The integral estimate is %4.3f.\n',I_new1)
    I_comp1 = (16*I_new1 - I_old1)/15; % Richardson's extrapolation from two estimate with an error order of h^4
    error1 = abs(I_comp1 - I_old1); % calculate the error
    fprintf('The estimate of the error is %4.3f.\n', error1)
    if error1 <= 0.005 % if the error is within the absolute error of 0.005, break the while loop
        break
    end
    I_old1 = I_new1; % reassign I_old for next evaluation
    h1 = h1/2; % decrase the size of h
end
fprintf('After evaluate bat_function_1 %d times, the final intergral for funtion_1 is %4.5f.\n\n', times1, I_new1)

% calculate the second function
I_old2 = 0;
times2 = 0;
h2 = h;
fprintf('Evaluating bat_function_2.\n')
while 1
    times2 = times2 + 1; % use times to record the number of times to evaluate
    fprintf('Attempt %d.\n', times2)
    I_new2 = run_simp(@bat_func_2, a, b, h2); % call simpson's 1/3 method
    fprintf('The integral estimate is %4.3f.\n',I_new2)
    I_comp2 = (16*I_new2 - I_old2)/15; % Richardson's extrapolation from two estimate with an error order of h^4
    error2 = abs(I_comp2 - I_old2); % calculate the error
    fprintf('The estimate of the error is %4.3f.\n', error2)
    if error2 <= 0.005 % if the error is within the absolute error of 0.005, break the while loop
        break
    end
    I_old2 = I_new2; % reassign I_old for next evaluation
    h2 = h2/2; % decrase the size of h
end
fprintf('After evaluate bat_function_2 %d times, the final intergral for funtion_2 is %4.5f.\n\n', times2, I_new2)

% calculate the total area of the bat sign
area = I_new1 - I_new2;
fprintf('The total area of the bat sign is %4.3f.\n', area)

%% Problem 2
a2 = 3;
b2 = 7;
end

function [I] = run_simp(f, a, b, h) % simpson' 1/3 method
N = (b-a)/h; % calculate the number of points
fprintf('The number of points is %d.\n', N);
xi = a:h:b; % create xi for evaluate
I = h/3*(f(xi(1))+2*sum(f(xi(3:2:end-2)))+4*sum(f(xi(2:2:end)))+f(xi(end))); % 1/3 method
end

function val = bat_func_1(xv) % copy bat__func_1 function
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
end 

function val = bat_func_2(xv) % copy bat__func_2 function
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
end