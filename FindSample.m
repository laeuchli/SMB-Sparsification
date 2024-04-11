function [y] = FindSample(p,q,n,startE,stopE,stepE, scale,nnM)

Evals = [];
Evals2= [];
percent =[];
totality = [];
for(E=startE:stepE:stopE)
    
    
    denom = 2*q*(n/2)^2+2*p*(n/2)^2;
    pp= (p*(n/2)^2*2)/denom;
    pq= (q*(n/2)^2*2)/denom;
    pp= (pp/(pp+scale*pq));
    pq=1-pp;
    pp=(E*pp)/((n/2)^2*2);
    pq=(E*pq)/((n/2)^2*2);
   
    a=pp*n;
    b=pq*n;

   % Evals2= [Evals2 (pp-pq)/sqrt(pp)-sqrt(n*log(n))/(pp*((n/2)^2*2))];
    Evals = [Evals (a-b)^2-2*(a+b)]
    totality = [totality E/nnM];
    alpha2= ((a-b)^2-2*(a+b))/((a-b)^2);
    if(alpha2> 0)
         percent= [percent 1/2*(1+erf(sqrt(alpha2/(2*(1-alpha2)))))]
    else
        percent = [percent 0 ];
    end
end

y = Evals;

end