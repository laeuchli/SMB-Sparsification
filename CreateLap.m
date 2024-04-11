function [y] = CreateLap(A)
newA=-A;

for(i=1:size(A,1))
    
    newA(i,i)=0;
    
end

for(i=1:size(A,1))

    newA(i,i)=-sum(newA(:,i));

end

y = newA;
