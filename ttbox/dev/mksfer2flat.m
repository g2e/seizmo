function [zf,varargout]=mksfer2flat(rp,zs,varargin)
%mksfer2flat    flattening earth transform: spherical to flat earth
%
% call: [zf,varargout]=mksfer2flat(rp,zs,varargin);
%
%               rp: planetary radius
%               zs: depths below surface in spherical planet
%         varargin: variable number of inputs (eg. vp,vs,etc) that are a
%                   function of spherical depth
%
%               zs and rp must be in the same units.
%               varargin inputs need to be the same length of vs
%
% result:           zf: equivalent depths in flattened planet
%            varargout: variable number of outputs corresponding to 
%                       varargin inputs transformed to work with flat earth
%                       depths.
%                       
%         Units are preserved. Planet center will be depth Inf.
%
% Traveltimes through the flattened model will be equal to the traveltimes
% in the spherical model.
%
% References: Mueller, G (1991), Inversionstheorie, Univ. Frankfurt
%
% Martin Knapmeyer, 18.04.2002, 10.11.2003, 10.02.2004
% Garrett Euler, 18.02.2008

%%% depth to inverse normalized radius
r=rp./(rp-zs);

%%% transformation of depths
zf=real(rp.*log(r)); % assure depths are real

%%% transformation of velocities or whatever
for i=1:nargin-2
    varargout{i}=varargin{i}.*r(:,ones(size(varargin{i},2),1));
end

end
