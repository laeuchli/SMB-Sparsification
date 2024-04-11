function [y1 y2  py1 py2] = TestData(L,DLD,rDLD,rL,minv,maxv,stepsize,minv2,maxv2,stepsize2, solver)

datagood =[];
databad =[];
sz= max(size(DLD));

p1 = [];
p2 = [];
nn=size(L,1);

for(i=minv:stepsize:maxv)
    sG=CreateSparseGraph(L, i, rDLD);
    A = sG;
    sG=CreateLap(real(logical(sG)));

    D=diag(diag(sG));
    
    dD=find(diag(D)==0);
   % A(dD,dD+1)=1;A(dD+1,dD)=1;
%     for(q=1:size(dD,1))
%         A(dD(q),mod(dD(q),10000)+1 )=1;
%         A(mod(dD(q),10000)+1,dD(q))=1;
%     end
    A=real(logical(A));
    sG=CreateLap(real(logical(A)));
    D=diag(diag(sG));   
    
    r= sqrt(sum(diag(D))/nn);
    O=ones(size(sG))*(r^2-1);
    avgD=sum(diag(D))/nn;
    %u1=ones(10000,1)*sqrt((r^2-1));%
    %tic;vq=powerfast2(r*sG,10^-8,25,u1);toc%
% %     tic;vs=powerfast(sG,10^-8,20);toc
% %     max( size(intersect(bb(1:5000),1:5000),1) ,size(intersect(bb(5001:10000),1:5000),1))
% %     [aa bb]=sort(vs(:,1));
  %  norm(D-avgD*speye(size(sG)))
 %   max(diag(D-avgD*speye(size(sG))));
 %   mean(abs(diag(D-avgD*speye(size(sG)))));
    %defL=O-r*A+D;
    tau = r;%5%(r^2-1);
    Or=(ones(size(sG))/nn)*tau;
    
    d=(diag(sG))+tau;;
    d=d.^(-1/2);
    DD=diag(d);
    D=diag(diag(D)+tau);
   % reg=DD*(Or+A-D)*DD;
    regA=Or+A-D;
   % [v3 e3]=eigs(sG,2,'sm');
    [v3 e3]=eigs(regA,D,2,'sm');
  % [v3 e3]= eig(regA, D);
% % %     
% % %     %[v3a e3a]=eigs(diag(diag(D)+tau)^-1*regA+ortho*ortho',2,'sm');
% % %     [v3 e3]=eigs(regA,diag(diag(D)+tau),2,'sm');
% % %     [v e]=eig(reg);
% % %     v1=v(:,2);
% % % 
% % %     D=diag(diag(sG));
% % %     refDLD=DD*(A-D)*DD;
% % %     refDL= D^-1*(A-D);
% % %     refL=A-D;
% % % 
% % %     
% % %     d=(diag(sG));;
% % %     d=d.^(-1/2);
% % %     DD=diag(d);
% % %     [v e]=eig(full(refDLD));
% % %     [v3 e3]=eigs(refL,D,2,'sm');
% % %     v1=v(:,9999);
% % %     vv3=v3(:,2);
% % %     norm(D.^(1/2)*vv3-v1);
% % %    defL=O-r*(A-D)+(1-r)*D;
% % %    % defL=O-r*(A-D)+(1-r)*avgD*speye(size(sG)) ; 
% % %     defLalt=O-r*(A-D) ; 
% % %     mean(abs(diag(D-avgD*speye(size(sG)))))
% % %      mean(abs(diag(D-avgD*speye(size(sG)))))/mean(diag(D))
  % [tv te]=eig(full(defL));
    if(solver==1) 
        
        u1=ones(10000,1)*sqrt((r^2-1));%
       % tic;vq=powerfast2(r*sG,10^-8,25,u1);toc%
   % tic;[t2v te2]=eigs(defLalt,1,'sa');;toc
   zz2=ones(10000,1);
   zz2=zz2/norm(zz2);
%   tic; t2v=powerfast4(A-D,D,100,ones(10000,1)*sqrt(tau/nn),zz2); toc
  t2v=v3(:,2);
   % t2v=vq;
     nnz(sG)
   % sG=CreateLap(real(logical(sG-diag(diag(sG)))));
    %[v e]=eig(full(sG));

    %[v e]=eigs(sG,2,'sm');
   % [v e]=eigs(sG+1/sqrt(sz)*ones(sz,1)*ones(sz,1)',1,'sm');
   % [aa bb]=sort(v(:,2));
    [aa bb]=sort(t2v(:,1));
     datagood = [datagood max( size(intersect(bb(1:5000),1:5000),1) ,size(intersect(bb(5001:10000),1:5000),1))]
     
    p1= [p1 nnz(sG)/39997804]
    else
         [v e]=eigs(sG,2,'sm');
         [aa bb]=sort(v(:,2));
         datagood = [datagood max( size(intersect(bb(1:5000),1:5000),1) ,size(intersect(bb(5001:10000),1:5000),1))]
    end
end



for(i=minv2:stepsize2:maxv2)
   
    sG=CreateSparseGraph(L, i, rL);
    A = sG;
    sG=CreateLap(real(logical(sG)));

    D=diag(diag(sG));

    
    dD=find(diag(D)==0);
   % A(dD,dD+1)=1;A(dD+1,dD)=1;
    for(q=1:size(dD,1))
        A(dD(q),dD(q)+1)=1;
        A(dD(q)+1,dD(q))=1;
    end
    A=real(logical(A));
    sG=CreateLap(real(logical(A)));
    
%     r= sqrt(sum(diag(D))/nn);
%     avgD=sum(diag(D))/nn;
%     O=ones(size(sG))*(r^2-1);
%     defL=O-r*A+D;
%     
%     
     d=(diag(sG));;
     d=d.^(-1/2);
     DD=diag(d);
     DsGD=DD*sG*DD;
   % [tv te]=eig(full(defL));
     if(solver==1) 
    [t2v te2]=eigs(DsGD,2,'sa');
     nnz(sG)
  %  [v e]=eig(full(sG));
%    [v e]=eigs(sG+1/sqrt(sz)*ones(sz,1)*ones(sz,1)',1,'sm');
  %  [v e]=eigs(sG,2,'sm');
   % [aa bb]=sort(v(:,2));
    [aa bb]=sort(t2v(:,2));
    databad = [databad max( size(intersect(bb(1:5000),1:5000),1) ,size(intersect(bb(5001:10000),1:5000),1))]
     p2= [p2 nnz(sG)/39997804]
     else
         [v e]=eigs(sG,2,'sm');
         [aa bb]=sort(v(:,2));
          databad = [databad max( size(intersect(bb(1:5000),1:5000),1) ,size(intersect(bb(5001:10000),1:5000),1))]
         
     end
end





datagood



y1 = datagood;
y2 = databad;
py1 = p1;
py2 = p2;

end