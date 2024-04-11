function [y] = powerfast2(A,tol,iter,u)

inputType = 'laplacian';        % The input matrix A is a graph Laplacian
solver = 'lamg';                % Or 'cmg', or 'direct'
lamg    = Solvers.newSolver(solver, 'randomSeed', 1);
setup   = lamg.setup(inputType, A);
b = rand(size(A,1), 1);
b = b - mean(b);    
for(i=1:iter)
 
    x=shermansolver(A,b,u,lamg); % [x, ~, ~, details] = lamg.solve(setup, b, 'errorReductionTol', 1e-8);
    b=x/norm(x);
end

y = b;

end