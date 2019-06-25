% MCEN 3030
% PROJECT 5
% HANWEN ZHAO
% MEID: 650-703
function [] = CP5_650703() % main function
% Problem 1
a = 0; b = 4; h = 0.1; yINI = 0;
[t1,y1] = Euler_Explicit(@ProblemOneODE,a,b,h,yINI);
%semilogx(t1,y1)
plot(t1,y1)
%[t2,y2] = Euler_Implicit(@ProblemOneODE,a,b,h,yINI);
%semilogx(t2,y2)
%plot(t2,y2)
%[t,y] = ode45(@ProblemOneODE,a:h:b,yINI);
%plot(t,y)
end

function [dydx] = ProblemOneODE(x,y)
%dydx = -1.2*y + 7*exp(-0.3*x);
%dydx = 0.8*y^(3/2) + 10*2000*(1-exp(-3*x));
dydx = -200000*y + 199000*y^(2/3)*exp(-x) + exp(-x);
end

function [x,y] = Euler_Explicit(ODE,a,b,h,yINI) % Euler Forward
x(1) = a; y(1) = yINI;
N = (b-a)/h;
for i = 1:N
    x(i+1) = x(i) + h;
    y(i+1) = y(i) + ODE(x(i),y(i))*h;
end
end

function [x,y] = Euler_Implicit(ODE,a,b,h,yINI) % Euler Backward
x(1) = a; y(1) = yINI;
N = (b-a)/h;
for i = 1:N
    x(i+1) = x(i) + h;
    ynew = y(i) + h*(ODE(x(i),y(i)));
    y(i+1) = y(i) + ODE(x(i+1),ynew);
end
end

function [x,y] = odeRK4(ODE,a,b,h,yINI)
x(1) = a; y(1) = yINI;
N = (b-a)/h;
for i = 1:N
    x(i+1) = x(i)+h;
    K1 = ODE(x(i),y(i));
    K2 = ODE(x(i)+h/2,y(i)+K1*h/2);
    K3 = ODE(x(i)+h/2,y(i)+K2*h/2);
    K4 = ODE(x(i+1),y(i)+K3*h);
    y(i+1) = y(i) + (K1+2*K2+2*K3+K4)*h/6;
end
end

function [x,y] = odeRK5(ODE,a,b,h,yINI)
x(1) = a; y(1) = yINI;
N = (b-a)/h;
for i = 1:N
    x(i+1) = x(i)+h;
    K1 = h*ODE(x(i),y(i));
    K2 = h*ODE(x(i)+h/4,y(i)+K1/4);
    K3 = h*ODE(x(i)+h/4,y(i)+(K1+K2)/8);
    K4 = h*ODE(x(i)+h/2,y(i)+(-K2+2*K3)/2);
    K5 = h*ODE(x(i)+3*h/4,y(i)+(-3*K1+9*K4)/16);
    K5 = h*ODE(x(i)+h,y(i)+(-3*K1+2*K2+12*K3-12*K4+8*K5)/7);
    y(i+1) = y(i) + (7*K1+32*K3+12*K4+32*K5+7*K6)/90;
end
end


