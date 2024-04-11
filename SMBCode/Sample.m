function [y] = Sample(p,  n )

    %[vals perm]=sort(p);
    %iperm(perm)=1:max(size(p));
    p=(p-min(p))/sum((p-min(p)));
   ip = p;
    % ip=max(p)-p;
    y= randsample(1:max(size(p)),n,true,ip);
  %  randvals =  randpdf(1-vals,1:max(size(p)),[n 1]);
   % randvals = ceil(randvals);
   % randvals(find(isnan(randvals)))=[];
   % y = iperm(randvals);

end