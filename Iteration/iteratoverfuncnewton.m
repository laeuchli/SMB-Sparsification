function [y y2] = iteratoverfuncnewton()
zz=-10:.01:10;
zz(1001)=[];
allvals = [];
correctedvs =[];
c=2
for(i=1:max(size(zz)))
   
    values = newton(zz(i),100);
    correctedvs = [correctedvs -(imag(c*values^2/(pi*zz(i)))) ];
    allvals =[ allvals values];
    
end

y= allvals;
y2 = correctedvs;