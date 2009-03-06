function [vs,zs]=mkflat2sfer(vf,zf,rp);
% mkflat2sfer.........flattening earth transform: flat to spherical earth
%
% call: [vs,zs]=mksfer2flat(vf,zf,rp);
%
%               vf: velocities in flat Planet at depths given in zf
%               zf: depths below surface in flat planet, at which velocites
%                   are given
%               rp: planetary Radius
%
%               zf and rp have to be in the same units.
%               vf and zf may be matrices of arbitrary size.
%
% result: vs: equivalent velocity in spherical Planet
%         zs: equivalent depths in spherical Planet
%
%         vs and zs will be in the same units as vf and zf are.
%         Depth inf will be transformed into the center of the Planet.
%
% wave traveltimes through the spherical model will be equal to the traveltimes
% in the flat model.
% References: Mueller, G (1991), Inversionstheorie, Univ. Frankfurt
%
% Martin Knapmeyer, 18.04.2002, 10.02.2004


%%% transformation of depths into distance from center
r=exp(-zf/rp); %r=rp*exp(-zf/rp);


%%% transformation of velocity
vs=r.*vf; %vs=r.*vf/rp;

%%% transformation of radii into depths beneath surface
zs=rp*(1-r); %zs=rp-r;


%%% that's all, folks.
 


