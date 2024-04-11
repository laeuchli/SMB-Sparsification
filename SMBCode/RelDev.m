function [y]= RelDev(pInv,L, i , j)

reldev = [];
for( ii=1:max(size(i)))

        R=effrefij(pInv,i(ii),j(ii));
        relR = (R-1/L(i(ii),i(ii))-1/L(j(ii),j(ii)));
        reldev = [reldev  relR];

end

y = reldev;
