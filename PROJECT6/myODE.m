function [ phi ] = myODE( t,y )

m=1.5;  %set constants
g=9.81;
c=0.5;

% dy(1)/dt = phi1 = y(2) = velocity
% dy(2)/dt = phi2 = -g +c/m*y(2)^2
% phi = [phi1 ; phi2] only need to input RHS (phi) for each ODE where phi1
% corresponds to the derivative of y(1) and phi2 with dy2/dt

if y(2) > 0  %drag is in negative y when moving up
    phi = [y(2); -g - c/m*(y(2))^2];
else         %other wise it is positive (none when v=0)
    phi = [y(2); -g + c/m*(y(2))^2];
end

end

