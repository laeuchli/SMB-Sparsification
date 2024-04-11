function [y] = func1examplenewt(z, guess)
c=2;
y=guess*z-c*guess^2-1;
%y=guess-(1/(z-c*guess));