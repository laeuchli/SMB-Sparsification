function [y] =newton(z,iter)
c=2;

guess =  -i;%(z-sqrt(z^2-4*c))/(2*c);%z*i+z;
current = 0;
for(ii=1:iter)
    f= func1examplenewt(z,guess);
    fp = z-2*c*guess;
    %fp = (-guess)/(z-2*c*guess);
  %  fp = 1/8*(z/(sqrt(z^2-8))+1);
    current = guess - f/fp;
     
    guess = current;
    
end
y = current;
