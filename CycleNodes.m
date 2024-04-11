function [y]= CycleNodes(L,k)

data=[];
for(i=1:max(size(L)))
    data=[data CycleEst(L,k,i)];
end
y= data;
end