function [y] = shermansolver2(A,b,u,lamg,setup,lowerterm)


[x, ~, ~, details] = lamg.solve(setup, b, 'errorReductionTol', 1e-8);
%x=A\b;%[x, ~, ~, details] = lamg.solve(setup, b, 'errorReductionTol', 1e-8);
tempb=dot(u,x);
[x2, ~, ~, details] = lamg.solve(setup, u*tempb, 'errorReductionTol', 1e-8);
%x2=A\(u*tempb);
y=x-(x2)/(lowerterm+1);

