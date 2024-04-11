%k1=[ones(2500,1);ones(2500,1)*.6];
k1=ones(5000,1);
k2=ones(5000,1);
k1=k1*200;
k2=k2*120;
theta1 = ones(5000,1)*50;
theta2= -ones(5000,1)*50;

k=[k1 theta1;k2 theta2];

m = norm(sum(k))/2;


Ab = Abracket(k,m);
A = poissrnd(Ab);% Ainit(Ab,k);

L=CreateLap(A);

d=(diag(L));
d=d.^(-1/2);
D=diag(d);
DLD=D*L*D;
disp('Model generated')
pDLD=pinv(DLD);
pL=pinv(L);
disp('Inverse generated')