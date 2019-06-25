clear all

y_init = 0; %IC
h = 0.1;
t_end = 4;
t = 0:h:t_end;
imax = length(t);

%RHS of eqn being solved, phi
f = @(t,y) -200000*y + 199000*y^(2/3)*exp(-t) + exp(-t);

%Let x replace y (x=y) when iterating through true y value in Newton's method
%once we converge to a given yi, we'll move to the next y value
%Bring y0 to the other side of eqn to match g(x) = 0 form for Newton's
g = @(x,y,t,h) x + 200000*x*h - 199000*x.^(2/3)*exp(-t)*h - exp(-t)*h - y;
%Take derivative of this g(x) with respect to x 
g_p = @(x,y,t,h) 1 + 200000*h - 199000*(2/3)*x^(-1/3)*exp(-t)*h;

% Euler (will diverge for this problem)
% y = zeros(1,imax);
% y(1) = y_init;
% for i = 1:imax-1
%     y(i+1) = y(i) + f(t(i),y(i))*h;
% end


% Implicit (Backward) Euler
y_imp = zeros(1,imax);
y_imp(1) = 0.0; %need reasonable initial condition for Newton's method
for i = 1:imax-1
    if i == 1, x=0.7; else  x = y_imp(i); end %can't divide by zero 
    %so don't start Newton iteration at y0 to determine y1, instead pick
    %another guess
    while true
        x_new = x - g(x,y_imp(i),t(i+1),h) / g_p(x,y_imp(i),t(i+1),h);
        if abs((x_new-x)/x) < 100*eps  %once we converge
            break
        else
            x = x_new;
        end
    end
    y_imp(i+1) = x;     %set our y value
end


%--------------------------------------------------------------------------
figure(1)
semilogx(t,y_imp,'r.-') %must use implicit Euler to ensure convergence
%axis([1E-4 1E4 -1 1])
figure(2)
plot(t,y_imp,'r.-') %must use implicit Euler to ensure convergence