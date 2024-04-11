function [a,b]=EstAB(A,K)
    n=max(size(A));
    d=2*nnz(A)/n;
    kn=3;%floor(log(n)^(1/4));
    f=(2*kn*K-d^kn)^(1/kn);
    a=d+f;
    b=d-f;

end