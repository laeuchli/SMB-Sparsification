function [R]= ChangeRands(e,v,range1,range2,ee1,ee2)

sumnoeigRDLD = 0;
for(i=range1:range2 )

    sumnoeigRDLD=sumnoeigRDLD + 1/e(i,i)*(ee1'*v(:,i)*v(:,i)'*ee2);
    
end

R=sumnoeigRDLD;

end