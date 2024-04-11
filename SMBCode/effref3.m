function [y]= effref3(pInv, i , j)

er =zeros(max(size(i)),1);
for( ii=1:max(size(i)))


        er(ii)=effrefij(pInv,i(ii),j(ii));  %=[er effrefij(pInv,i(ii),j(ii))];

end

y = er;
