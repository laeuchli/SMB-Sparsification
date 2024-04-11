n=1000;
AA=sparse(n,n);
p=[.02 .005;.005 .02];
for(i=1:n)
    for(j=1:n)
        
        if(rand < AZ(i,j))
            AA(i,j)=1;
        end
    end
end

c1=zeros(10000,1);
c1(1:5000)=1;
p=[.2 .05;.05 .2];
A=generateSbm(c1+1,p);
L=CreateLap(A);

d=(diag(L));
d=d.^(-1/2);
D=diag(d);
DLD=D*L*D;
disp('Model generated')
pDLD=pinv(DLD);
pL=pinv(L);
disp('Inverse generated')


[i j z]=find(A);
rDLD=effref3(pDLD,i,j);
rL=effref3(pL,i,j);
disp('Resistance generated')
[v1 e1]=eig(A);
[v2 e2]=eig(L);
[v3 e3]=eig(DLD);

%fun2=@(x) x.*((n./((cin+cout)*pi))*sqrt(2*(cin+cout)-x.^2))
%(integral(fun2,0,19.99)*2-sum(abs(eee1)))/sum(abs(eee1))
[vz ez]=eig(DLD-eye(size(L)));
eez=diag(ez);
eeez=eez(3:10000);
er=max(eeez)
n=10000
eee3=diag(e3);
eee1=1./eee3(3:10000);
 %(evals(1)/1.75)^102
 %evals3=(evals)/1.7518;
%evals=1./((eee3(3:10000)-.29));
a=.5*n;
b=.1*n;
%1-sqrt(2*(a + b))/dd
%1+sqrt(2*(a + b))/dd
%1+e1(9998,9998)/dd
%1-e1(9998,9998)/dd
er=((1+e1(9998,9998)/dd)- (1-e1(9998,9998)/dd)) /2;
fun2=@(x) x.*((2*n./(er^2*pi))*sqrt(er^2-x.^2))
(integral(fun2,0,er)*2-sum(abs(eee1)))/sum(abs(eee1))

 fun2=@(x) (x+.9487+er).*((2*n./(er^2*pi))*sqrt(er^2-x.^2))
(integral(fun2,-er,er)-sum(abs(eee1+1)))/sum(abs(eee1+1))

%  fun2=@(x) (x+0.9761+er).*((2*n./(er^2*pi))*sqrt(er^2-x.^2)) %shifted
%  specturm
% (integral(fun2,-er,er)-sum(abs(eee1+1)))/sum(abs(eee1+1))

fun2=@(x) (1./(x+.9487+er)).*((2*n./(er^2*pi))*sqrt(er^2-x.^2))
(integral(fun2,-er,er)-sum(abs(eee1)))/sum(abs(eee1))

pows=12;
fun2=@(x) ((1./(x+.9728+er)).^pows).*((2*n./(er^2*pi))*sqrt(er^2-x.^2))
(integral(fun2,-er,er)-sum(abs(eee1.^pows)))/sum(abs(eee1.^pows))
%fun2=@(x) (((1./(x+.9761+er-.29))/1.7518).^pows).*((2*n./(er^2*pi))*sqrt(er^2-x.^2))
fun2=@(x) (((1./(x+.9728+er))/ 3.0039).^pows).*((2*n./(er^2*pi))*sqrt(er^2-x.^2))

sumnoeigR=[];
sumnoeigS=[];
sumgR=[];
sumgS=[];
%
sumnoeigRL=[];
sumnoeigSL=[];
sumgRL=[];
sumgSL=[];
%
sumnoeigRDLD=[];
sumnoeigSDLD=[];
sumgRDLD=[];
sumgSDLD=[];
trueR=[];

for(k=1:100)
    ee1=zeros(10000,1);
    ee2=zeros(10000,1);
    k
    ii=randi(500)+5000;
    jj= randi(500)+50
    ee1(ii)=1;
    ee2(jj)=1;
    if( A(ii,jj)==1)
        trueR=[trueR ee1'*pDLD*ee2];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sumnoeigR = [sumnoeigR ChangeRands(speye(size(A)),v1,1,max(size(A))-2,ee1,ee2 )];
    
    sumnoeigS = [sumnoeigS ChangeRands(speye(size(A)),v1,999,999,ee1,ee2 )];
    
    
    sumgR = [sumgR ChangeRands(e1,v1,1,max(size(A))-2,ee1,ee2 )];
    
    sumgS = [sumgS ChangeRands(e1,v1,999,999,ee1,ee2 )];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sumnoeigRL = [sumnoeigRL ChangeRands(speye(size(A)),v2,3,max(size(A)),ee1,ee2 )];
    
    sumnoeigSL = [sumnoeigSL ChangeRands(speye(size(A)),v2,2,2,ee1,ee2 )];
    
    
    sumgRL = [sumgRL ChangeRands(e2,v2,3,max(size(A)),ee1,ee2 )];
    
    sumgSL = [sumgSL ChangeRands(e2,v2,2,2,ee1,ee2 )];
    
    %ee2'*pinv(L)*ee1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sumnoeigRDLD=[sumnoeigRDLD ChangeRands(speye(size(A)),v3,3,max(size(A)),ee1,ee2 )];
    
    sumnoeigSDLD = [sumnoeigSDLD ChangeRands(speye(size(A)),v3,2,2,ee1,ee2 )];
    
    
    sumgRDLD = [sumgRDLD ChangeRands(e3,v3,3,max(size(A)),ee1,ee2 )];
    
    sumgSDLD = [sumgSDLD ChangeRands(e3,v3,2,2,ee1,ee2 )];
    
end
mean(sumnoeigR)
mean(sumnoeigS)
mean(sumgR)
mean(sumgS)
%
mean(sumnoeigRL)
mean(sumnoeigSL)
mean(sumgRL)
mean(sumgSL)
%
mean(sumnoeigRDLD)
mean(sumnoeigSDLD)
mean(sumgRDLD)
mean(sumgSDLD)
%ee=diag(e3);
 %2+1/n*1/ee(2)-2*(0+mean(sumgSDLD)) %  1.9998
%mean(zr(find(zr>0))) 1.9998
% ssr =[];
% for(k=1:1000)
% 
%     
%     ee1=randi(50);
%     ee2=randi(50);
%     
%      ssr = [ssr  RR(i(k),j(k))-1/d(i(k))-1/d(j(k))];
% 
% end

% for(si=1:max(size(i)))
%     
%     if(i(si) <= 50 && j(si) <=50)
%      ssr = [ssr  RR(i(si),j(si))-1/d(i(si))-1/d(j(si))];
%     end
%     
% end
%ee2'*pinv(DLD)*ee1