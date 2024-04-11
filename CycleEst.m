function [y]=CycleEst(L,k,v)

Lp=L-diag(diag(L));
A=real(logical(L-diag(diag(L))));

e=sparse(max(size(A)),1);
e(v)=1;

for (i=1:k)
    e=A*e;
end

D=abs(sum(Lp(find(e),find(e))));
D=diag(D);
newA=A(find(e),find(e));

A0=speye(size(newA));
A1=newA;


ne=sparse(max(size(newA)),1);
ne(find(e)==v)=1;

for(i=1:k)
%     te1=newA*ne;
%     te1=A1*ne;
%     te2=A0*ne;
%     te2=D*te2;
%     ne=te1-te2;
tempA=A1*newA-D*A0;
A0=A1;
A1=tempA;
  
end

y = A1(find(e)==v,find(e)==v);