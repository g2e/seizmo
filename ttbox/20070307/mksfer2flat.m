function [vf,zf]=mksfer2flat(vs,zs,rp);
% mksfer2flat.........flattening earth transform: spherical to flat earth
%
% call: [vf,zf]=mksfer2flat(vs,zs,rp);
%
%               vs: velocities in spherical Planet at depths given in zs
%               zs: depths below surface in spherical planet, at which velocites
%                   are given
%               rp: planetary Radius
%
%               zs and rp have to be in the same units.
%               vs and zs may be matrices of arbitrary size.
%
% result: vf: equivalent velocity in flattened Planet
%         zf: equivalent depths in flattened Planet
%
%         vf and zf will be in the same units as vs and zs are.
%         The center of the planet will be transformed to depth inf.
%
% wave traveltimes through the flattened model will be equal to the traveltimes
% in the spherical model.
% References: Mueller, G (1991), Inversionstheorie, Univ. Frankfurt
%
% Martin Knapmeyer, 18.04.2002, 10.11.2003, 10.02.2004

%%% transform spherical depth into Radius
r=rp-zs; 

%%% find zero elements of r = these cannot be mapped using the log
%%% r=0 will be mapped to zf=inf explicitly.
indies=find(r==0);
r(indies)=-1; % r=-1 is physically impossible and produces no warning

%%% transformation of dephts
zf=rp.*log(rp./r);
zf(indies)=abs(zf(indies))*inf;  % repairs the r=-1 substitution from above
zf=real(zf); % even flat depths are real numbers. MK10112003

%%% transformation of velocities
vf=vs./r;
vf=vf*rp;
vf(indies)=vf(indies)*inf; % repairs the r=-1 substitution from above

%%% that's all, folks. 


