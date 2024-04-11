function [y]= effref2(pInv)
n = max(size(pInv));

ertotal = [];
    for(i=1:n)
        er=[];
        for(j=1:n)
            
            if(i==j)
                continue;
            end
            e1 = zeros(n,1);
            e2 = zeros(n,1);
            e1(i)=1;
            e2(j)=1;
            ef = (e1-e2)'*pInv*(e1-e2);
            er=[er ef];
        end
        ertotal=[ertotal; er];
    end


y = ertotal;


end