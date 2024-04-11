function [y] = shermansolver(A,b,u,lamg,setup)



[x, ~, ~, details] = lamg.solve(setup, b, 'errorReductionTol', 1e-8);
tempb=dot(u,x);
[x2, ~, ~, details] = lamg.solve(setup, u*tempb, 'errorReductionTol', 1e-8);

y=x+x2;

