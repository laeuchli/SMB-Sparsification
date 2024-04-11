% Wigner semicirle law : 
% Probability distribution of the eigenvalues of symmetric centered random matrices
% as the limit when the dimension N tends to Infinity.
% (c)KHMOU Youssef,  Random Matrix theory, May, 2014.

N=10000;
%r=sqrt(N);
M=randn(N);
M=(M'+M)/sqrt(2*N); % M^T=M.

vcs=zeros(N,1);
for(i=1:N)
    vcs(i)=exp(i/5000);
end

vcs = vcs/norm(vcs);

P= diag(vcs);
M=P*M*P;

L=eig(M);

fun2=@(x) x.*((N./(2*pi))*sqrt(4-x.^2))
integral(fun2,0,2)*2
sum(abs(L))

Nbins=60;
[Y,X]=hist(L,Nbins);
Y=Y/sum(Y);
figure; h=bar(X,Y);
set(h,'FaceColor','none');
hold on;

% Theoretical Probability density function
x=linspace(-max(L)*sqrt(N),max(L)*sqrt(N),10000);
r=max(L);
alpha=1.5; % Scalling parameter, it depends on variable Nbins.
f=(2*alpha/(pi*(r^2)))*sqrt((r^2)-((x.^2)));
plot(x,f,'r','LineWidth',2);
axis([min(L)-r max(L)+r 0 max(Y)+max(Y)/4]);
xlabel('\lambda','FontWeight','bold','FontSize',12)
ylabel('PDF(\lambda)','FontWeight','bold','FontSize',12);
title('PDF the  Eigenvalues of centered random symmetric matrix');
Info='$f(\lambda)={2\over{\pi r^{2}}} \sqrt{{r^{2}}-{{\lambda^{2}}}}$'; 
text('Position',[-2*sqrt(N) max(Y)+max(Y)/8],'Interpreter','latex','String',Info,'FontSize',16);2^9