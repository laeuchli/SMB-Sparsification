function [y] = func1example(z, guess)
c=2;
y= guess*z-c*guess^2-1+guess;
%y=(c*guess^2+1)/z;
%y=(1/(z-c*guess));