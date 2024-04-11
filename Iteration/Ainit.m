function [y] = Ainit(A,k)

trueA = zeros(size(A));
for( i = 1: max(size(A)))
    for(j=1:max(size(A)))
        trueA(i,j) = poissrnd(A(i,j));
    end
end

y = trueA;
end