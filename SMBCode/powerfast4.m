function [y] = powerfast4(A,B,iter,u,orthov)

inputType = 'sdd';%'laplacian';        % The input matrix A is a graph Laplacian
solver = 'lamg';                % Or 'cmg', or 'direct'
lamg    = Solvers.newSolver(solver, 'randomSeed', 1);
setup   = lamg.setup(inputType, A);
b = rand(size(A,1), 1);
b = b - mean(b);    
D=diag(B);
D=1./D;
D=diag(D);

[x, ~, ~, details] = lamg.solve(setup, u, 'errorReductionTol', 1e-8);
lowerterm= dot(u,x);
for(i=1:iter)
    b=b-orthov*dot(orthov,b);
    b=b/norm(b);
    b=B*b;
    x=A(A,b,u,lamg,setup,lowerterm); % [x, ~, ~, details] = lamg.solve(setup, b, 'errorReductionTol', 1e-8);
    if(i==95)
        i=i;
    end
%     x=x-orthol*orthov*dot(orthov,b);
    x=x-orthov*dot(orthov,x);
    b=x/norm(x);
    
end

y = b;

end