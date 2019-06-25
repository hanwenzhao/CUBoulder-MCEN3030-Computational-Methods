close all
clear all
%send in ODE, time range, ICs for each ODE
[t45, y45] = ode45(@myODE, [0 1.6], [1, 10]);
plot(t45, y45, '*')
legend('y', 'v')

%Look up 2 points that bracket y=0, near t=1.6
%Perform linear interpolation to get the ~time where y=0