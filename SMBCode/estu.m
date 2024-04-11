function [y] = estu(L,i,j)


eu = [];
for( ii=1:max(size(i)))

      u = 2*L(i(ii),j(ii))/(L(j(ii),j(ii))*L(i(ii),i(ii)))-L(i(ii),i(ii))/(L(i(ii),i(ii))^2)-L(j(ii),j(ii))/(L(j(ii),j(ii))^2);
      eu = [eu u];

end

y = eu;

end