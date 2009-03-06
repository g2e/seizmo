function [dny,dxn]=mkderive(x,y,n);
% mkderive.........n-th derivative of discrete function
%
% call: [dny,dxn]=mkderive(x,y,n);
%
%       x: x values of function samples
%       y: y values of function samples
%       n: the n-th derivattive will be computed
%
% result: dny: n-th differences of y
%         dxn: n-th differences of x
%
%         dny/dxn is an approximation to the n-th derivative of the function y(x)
%
% This routine approximates the derivative by iterative application of
% sample differences. Samples do _not_ need to be equidistant.
%
% Martin Knapmeyer, 01.06.2006


dny=y;
dxn=x;
for indy=1:n
    dny=diff(dny)./diff(dxn);
    dxn=dxn(1:(end-1));
end; % for indy



