function [cp,cg]=swell_velocity(l,d)
%SWELL_VELOCITY    Returns speed of a swell given wavelength & ocean depth
%
%    Usage:    [c,g]=swell_velocity(lambda,depth)
%
%    Description:
%     [C,G]=SWELL_VELOCITY(LAMBDA,DEPTH) returns the phase velocity C &
%     group velocity G given the wavelength LAMBDA and water depth DEPTH.
%     Note that inputs are in meters and outputs in meters/second.
%
%    Notes:
%     - Handles arbitrary water depth and wavelength (includes capillary
%       terms).
%
%    Examples:
%     % What is the period?
%     c=swell_velocity(l,1000);
%     T=l/c;
%
%    See also: SWELL_BACKPROJ, SWELL_FORWARD

%     Version History:
%        May  16, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  16, 2012 at 15:05 GMT

% todo:

% check number of inputs
error(nargchk(2,2,nargin));

% require equal size or scalar
if(~isscalar(l) && ~isscalar(d) && ~isequal(size(l),size(d)))
    error('seizmo:swell_velocity:badInput',...
        'LAMBDA & DEPTH must be equally sized or scalars!');
end

% calculate phase & group speed at an
% arbitrary depth & arbitrary wavelength
tau=2*pi;
k=tau./l;
g=9.81;
gamma=7.25e-2; % room temp
rho=1000;
g=g+gamma/rho*k.^2;
thkd=tanh(k.*d);
cp=sqrt(g./k.*thkd);
cg=.5*cp.*(1+k.*d.*((1-thkd.^2)./thkd));

end
