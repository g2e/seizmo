function rho=mkbirch(vp);
% mkbirch...............calculate density from P velocity unsing Birch's law
%
% call: rho=mkbirch(vp);
%
%           vp: vector or matrix of p wave velocities, in m/sec
%
% result: rho: the densities corresponding to the velocities in VP.
%              RHO and VP are of equal size.
%              RHO will be in kg/m^3
%
% Birch's Law is a approxmation of type vp=a+b*rho, with empirical parameters
% a and b.
% This routine uses
%                          a=-0.98 km/sec
%                          b= 2.76 km*cm^3/(g*sec)
%
% (Birch's solution 5 for all rock samples except those with mean atomic weigth
%  above 24)
%
% reference:
% Birch, Francis (1961):
% The Velocity Of Compressional Waves In rocks To 10 Kilobars, Part 2;
% JGR, vol. 66, No. 7, 2199-2224
%
% Martin Knapmeyer, 29.03.1996

a=-0.98*1000; % dimension: meters per sec
b=2.76; % dimension: m*m^3/(kg*sec)

rho=(vp-a)/b;
