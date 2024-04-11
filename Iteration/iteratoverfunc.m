function [y y2] = iteratoverfunc()
zz=-10:.01:10;
zz(1001)=[];
allvals = [];
correctedvs =[];
c=2
for(i=1:max(size(zz)))
   
    values = fixedpoint(zz(i),1000);
    if(imag(values)>.1)
     a=1;
    end
    correctedvs = [correctedvs -(imag(c*values^2/(pi*zz(i)))) ];
    allvals =[ allvals values];
    
end

y= allvals;
y2 = correctedvs;