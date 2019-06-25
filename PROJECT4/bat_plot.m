clear all
clc
close all


xlabel('x')
ylabel('y')
title('original bat')

figure(2)
x = [-7:0.01:7];
plot(x,bat_func_1(x),'-r')
hold on
plot(x,bat_func_2(x),'--b')
xlabel('x')
ylabel('y')
title('reconstructed bat')
legend('b_1(x)','b_2(x)')
axis([-7.25 7.25 -5 5]);

