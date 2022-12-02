%Script that calculates the cubic spline interpolation function for a
%user defined experimental data set of x and f(x) values, and outputs the
%interpolated value of the function at a user defined x value.
%Written by Jonathan De Sousa
%Date: 02-03-2019

%clear workspace and screen
clear
clc

%Prompt user for experimental values
x=input('Enter x values as a column vector: '); 
f=input('Enter f(x) values as a column vector: ');
xval=input('Enter x value: ');

%Calculate cubic spline interpolant function 
polycoeff=coeff(x,f); %calculate polynomial coefficients of interpolant
yInterpVal=cubicEval(x,polycoeff,xval); %obtain interpolated f(x) value

%Display interpolated f(x) value at xval
disp(['Interpolated f(',num2str(xval),') value: ',num2str(yInterpVal)])

%end of script






