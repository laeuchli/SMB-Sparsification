function [y]= effrefij(pInv,i,j)
n = max(size(pInv));


er=[];
    %for(i=1:n)
     %   for(j=1:n)
            
%             if(i==j)
%                 continue;
%             end
            e1 = zeros(n,1);
            e2 = zeros(n,1);
            
            e1(i)=1;
            e2(j)=1;
            
            et= pInv(:,i)-pInv(:,j);
            et= et(i)-et(j);
        %    ef = (e1-e2)'*pInv*(e1-e2);
           % er=[er et];
      %  end
    %end


y = et;


end