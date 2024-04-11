function [y] = CreateSparseGraph(L, n, p)

QZ=L-diag(diag(L));
[i j z ] = find(QZ);

samples= Sample(p,n);


pvals=(p-min(p))/sum((p-min(p)));
%ipvals=max(pvals)-pvals;
%y = sparse(i(samples), j(samples), ones(max(size(samples)),1)./(n*pvals(samples)), size(L,1), size(L,2))+sparse(i(samples), j(samples), ones(max(size(samples)),1)./(n*pvals(samples)), size(L,1), size(L,2))';
y = sparse(i(samples), j(samples), ones(max(size(samples)),1), size(L,1), size(L,2))+sparse(i(samples), j(samples), ones(max(size(samples)),1), size(L,1), size(L,2))';

