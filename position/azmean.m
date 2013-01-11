function [avg]=azmean(az,dim)
%AZMEAN    Returns the mean azimuth of a set of azimuths
%
%    Usage:    avg=azmean(az)
%              avg=azmean(az,dim)
%
%    Description:
%     AVG=AZMEAN(AZ) returns the average azimuth of the azimuths in AZ.  AZ
%     is expected in degrees.  Operates down the first non-singleton
%     dimension.
%
%     AVG=AZMEAN(AZ,DIM) takes the mean across along dimension DIM.
%
%    Notes:
%
%    Examples:
%     % Average North scattered azimuths:
%     azmean(rand(100,1)-.5)
%
%    See also: AZDIFF

%     Version History:
%        Oct. 10, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 10, 2012 at 19:00 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% default dimension
if(nargin==1 || isempty(dim))
    dim=find(size(az)~=1,1);
    if(isempty(dim)); dim=1; end
end

% azimuths to unit vector, mean, unit vector to azimuth
avg=180/pi*atan2(mean(sind(az),dim),mean(cosd(az),dim));

end
