function pf=mkpsfer2flat(ps,rp);
% mkpsfer2flat...........transformation of ray parameter from spherical to flat s/rad -> s/km
%
% call: pf=mkpsfer2flat(ps,rp);
%
%          ps: ray parameter in sec/rad as needed fopr spherical calculations
%          rp: planetary radius
%
% result: pf: ray parameter in sec/km as needed for FLAT EARTH calculations
%
% Since the earth's sureface is transformed into a straight line of infinite length,
% the ray parameter has to be converted into an appropriate unit.
% 
%
% Martin Knapmeyer, 26.04.2002


%%% the transformation
pf=ps/rp; % is the same as (ps*pi/180)/(2*pi*rp/360); which is what really happens
