function [y] = fixedpoint(z,iter)
c=2;
if(z^2<4*c)
    c = 2;
end
guess = 0;

guess = -i*z;%;((z-sqrt(z^2-4*c))/(2*c));%z*i+z;

debug = [];
err = 0;
for(ii=1:iter)
       newguess = func1example(z,guess);
       err = [err  ((newguess)-(guess))];
       if(((newguess)-(guess)) < 10^-8)
           y = (guess);
           return;
       end
       guess = newguess;
       debug = [debug guess];
end
    y = (guess);
    

 