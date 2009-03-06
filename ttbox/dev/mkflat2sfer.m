function [zs,varargout]=mkflat2sfer(rp,zf,varargin)
%mkflat2sfer    flattening earth transform: flat to spherical earth
%
% call: [zs,varargout]=mksfer2flat(rp,zf,varargin);
%
%               rp: planetary radius
%               zf: depths below surface in flat planet
%         varargin: variable number of inputs (eg. vp,vs,etc) that are a 
%                   function of flat earth depth.
%
%               zf and rp must be in the same units.
%               varargin inputs need to be the same length of zf.
%
% result:           zs: equivalent depths in spherical planet
%            varargout: variable number of outputs corresponding to 
%                       varargin inputs transformed to work with spherical
%                       depths.
%                       
%         Units are preserved.  Center of earth will be depth rp.
%
% Traveltimes through the spherical model will be equal to the traveltimes
% in the flat model.
%
% References: Mueller, G (1991), Inversionstheorie, Univ. Frankfurt
%
% Martin Knapmeyer, 18.04.2002, 10.02.2004
% Garrett Euler, 18.02.2008


%%% transformation of depths into distance from center
r=exp(-zf/rp);

%%% transformation of radii into depths beneath surface
zs=rp*(1-r);

%%% transformation of velocity
for i=1:nargin-2
    varargout{i}=varargin{i}.*r(:,ones(size(varargin{i},2),1));
end

end
