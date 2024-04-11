function [y] = powerfast(A,tol,iter,orthov)

inputType = 'laplacian';        % The input matrix A is a graph Laplacian
solver = 'lamg';                % Or 'cmg', or 'direct'
lamg    = Solvers.newSolver(solver, 'randomSeed', 1);
setup   = lamg.setup(inputType, A);
b = rand(size(A,1), 1);
b = b - mean(b);    
for(i=1:iter)
    
    b=b-orthov*dot(orthov,b);
    b=b/norm(b);
    
    [x, ~, ~, details] = lamg.solve(setup, b, 'errorReductionTol', 1e-8);
    x=x-orthov*dot(orthov,x);
    b=x/norm(x);
    
end

y = b;

end