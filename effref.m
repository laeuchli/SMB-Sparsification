function [y]= effref(pInv)
n = max(size(pInv));


er=[];
    for(i=1:n)
        i
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
    end


y = er;


end