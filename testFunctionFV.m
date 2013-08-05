function [f] = testFunctionFV(x,shape)
% f = testFunctionFV(x,shape)
%
% Integrates the "testFunction()" over cells with edge positions
% given by x. This routine calls testFunction and numerically
% integrates the function over each cell using Boole's rule.
%
% f will have size(x)-1 cells.
%
% plot(0.05:.01:1, testFunctionFV(0:.01:1,'pulse') )

% Boole's rule
xa = x(:,1:end-1);
xb = x(:,2:end);
f0 = testFunction(xa,shape);
f1 = testFunction((3*xa+xb)/4,shape);
f2 = testFunction((xa+xb)/2,shape);
f3 = testFunction((1*xa+3*xb)/4,shape);
f4 = testFunction(xb,shape);
f=( 32*(f1+f3) + ( 7*(f0+f4) + 12*f2 ) )/90;
