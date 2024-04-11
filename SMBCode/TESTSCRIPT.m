c1=zeros(10000,1);
c1(1:5000)=1;
 p=[.5 .3;.3 .5];
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
rDLDc=effref3(pDLD,i,j);
rL=effref3(pL,i,j);
disp('Resistance generated')

%%%%%%%
tau = 34;%(r^2-1);
Or=(ones(size(A))/size(A,1))*tau;

d=(diag(A))+tau;
d=d.^(-1/2);
DD=diag(d);
D=diag(diag(L)); 
Dt=diag(diag(D)+tau);
% reg=DD*(Or+A-D)*DD;
regA=Or+A-Dt;

DtADt=Dt^(-1/2)*regA*Dt^(-1/2);
[v3 e3]=eigs(regA,D,2,'sm');
%%%%%%%%%%%%%%%%%%%%%


stepsize = 20000
basestart =20000

stepsize2 = 20000
basestart2 =20000
%[gd bd] = TestData(L,DLD,rDLD,rL,basestart,basestart+stepsize*25,stepsize);abs(rDLD-2*sum(eaa)/9998)
[gd bd] = TestData(L,DLD,abs(rDLD-2),abs(rDLD-2),basestart,basestart+stepsize*25,stepsize,basestart2,basestart2+stepsize2*25,stepsize2,1);
gdn=(gd)/5000
bdn=(bd)/5000

figure(1)
plot([1:max(size(gdn))],gdn,'x')
figure(2)
plot([1:max(size(bdn))],bdn,'o')