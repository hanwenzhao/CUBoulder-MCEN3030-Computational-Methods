%MCEN 303
%Project 2
%MEID: 650-703
%% Main function
function [] = cp2_650703()
clear All
clc


%a
fprintf('Problem a \n') % solve problem a
n = 5;
m = 10^2;
kmax = 10;
mainFun(n,m,kmax)

%b
fprintf('Problem b \n') % solve problem b
n = 5;
m = 10^2;
kmax = 50;
mainFun(n,m,kmax)

%c
fprintf('Problem c \n') % solve problem c
n = 5;
m = 10^6;
kmax = 50;
mainFun(n,m,kmax)

%d
fprintf('Problem d \n') % solve problem d
n = 20;
m = 1;
kmax = 10;
mainFun(n,m,kmax)

%e
fprintf('Problem e \n') % solve problem e
n = 20;
m = 0;
kmax =50;
mainFun(n,m,kmax)

end
% Problem #1
function [A,b] = buildMatrix(n,m) % build matrix function
A = magic(n) + eye(n) * m; % build A
for i = 1:n % build b
    b(i) = 1/i;
end
b = b';
end
% Problem #2
function [x] = status(A)
x = -999; % return -999 if matrix is not well conditioned and diagonally dominant
syms diaDom % varible to record diagonally dominant or not
n = length(A); 
for i = 1:n
        B = abs(A(i,:)); % B is single row of A
        if abs(A(i,i)) > sum(B) - abs(A(i,i)) % determine A is diagonally dominant or not
            fprintf('Diagonally dominant. \n')
            diaDom = 1;
            break
        else
            diaDom = 0;
        end
end

if cond(A) < 10^16 % if condition number is less than O(10^15)
    fprintf('Well conditioned. \n')
end
if diaDom == 1 && cond(A) < 10^15 % if A is diagonally dominant and the condition number is less than O(10^15)
    x = 0; % set status varible to 0
    fprintf('The status varible is 0. \n')
end
if diaDom == 0
    fprintf('Not Diagonally dominant. \n')
    x = -1;
end
if cond(A) >= 10^16 % if the condition number is greater than O(10^15)
    fprintf('Ill-conditioned. \n') % matrix is ill conditioned, set x to -1
    x = -2;
end
end
% Problem #3
function [x] = solveJacobi(A,b,kmax) % use jocobi iterative method
eps = 10^-10; % relative error tolerance
converged = 0; % record the solution converged or not
n = length(b);
p = zeros(n,1); % assume all initial value to zero
for k = 1:kmax
    for i = 1:n
        x(i) = (b(i) - A(i,[1:i-1,i+1:n])*p([1:i-1,i+1:n]))/A(i,i); % estimate the solution
    end
    err = abs(norm(x'-p)/norm(p)); % calculate error
    fprintf('The relative error is %d \n',err)
    if err < eps % if error less than tolerance, stop iternation
        fprintf('The total estimated error is less than a relative error tolerance. \n')
        converged = 1; % the solution cinverge
        break
    end
    p = x'; % transpose the matrix
end
x = x';
if converged == 0 % if the solution did not converge, set x to -999
    fprintf('The system did not converge before %d iterations are performed. \n', kmax)
    x = -999;
end
end
% Problem #4
function [x] = mainFun(n,m,kmax) % main function
[A,b] = buildMatrix(n,m); % call function to build matrix
statusVarible = status(A); % calculate the status varible
if statusVarible == 0 % if matrix is diagonally dominant and the well condition
    x = solveJacobi(A,b,kmax); % solve x by using jacobi method
else
    x = -999;
end

end

