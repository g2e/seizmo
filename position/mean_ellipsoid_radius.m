function [r]=mean_ellipsoid_radius(method,ellipsoid)
%MEAN_ELLIPSOID_RADIUS    Returns the mean radius of an ellipsoid
%
%    Usage:    r=mean_ellipsoid_radius()
%              r=mean_ellipsoid_radius(method)
%              r=mean_ellipsoid_radius(method,[a f])
%
%    Description:
%     R=MEAN_ELLIPSOID_RADIUS() returns the mean radius of the WGS-84
%     reference ellipsoid based on the formula (2A+B)/3 where A is the
%     semimajor axis and B is the semiminor axis.
%
%     R=MEAN_ELLIPSOID_RADIUS(METHOD) specifies the method for finding the
%     mean radius of an ellipsoid.  Valid methods are 'MEAN', 'AUTHALIC',
%     and 'VOLUMETRIC'.  AUTHALIC gives the radius of a sphere with an area
%     equal to that of the ellipsoid.  VOLUMETRIC gives the radius of a
%     sphere with an equal volume to that of the ellipsoid.  The default
%     method is 'MEAN' which uses the formula above.
%   
%     R=MEAN_ELLIPSOID_RADIUS(METHOD,[A F]) allows specifying the ellipsoid
%     parameters A (equatorial radius in kilometers) and F (flattening).
%     This is compatible with Matlab's Mapping Toolbox function ALMANAC.
%
%    Notes:
%
%    Examples:
%     % Compare the 3 mean radii for the WGS-84 ellipsoid:
%     mean_ellipsoid_radius('mean')
%     mean_ellipsoid_radius('authalic')
%     mean_ellipsoid_radius('volumetric')
%
%    See also: GEOGRAPHICLAT2RADIUS

%     Version History:
%        Nov. 13, 2009 - initial version
%        Nov. 16, 2009 - divide by zero bug fixed in sphere 2 authalic
%        Feb. 11, 2011 - mass nargchk fix
%        Feb.  9, 2012 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(0,2,nargin));

% default ellipsoid - WGS-84 Reference Ellipsoid
if(nargin<2 || isempty(ellipsoid))
    % a=radius at equator (major axis)
    % f=flattening
    a=6378.137;
    f=1/298.257223563;
else
    % manually specify ellipsoid (will accept almanac output)
    if(isnumeric(ellipsoid) && numel(ellipsoid)==2 && ellipsoid(2)<1)
        a=ellipsoid(1);
        f=ellipsoid(2);
    else
        error('seizmo:mean_ellipsoid_radius:badEllipsoid',...
            ['Ellipsoid must a 2 element vector specifying:\n'...
            '[equatorial_km_radius flattening(<1)]']);
    end
end

% get semiminor axis
b=a*(1-f);

% default method
if(nargin==0 || isempty(method)); method='mean'; end

% check method
if(~ischar(method))
    error('seizmo:mean_ellipsoid_radius:badMethod',...
        'METHOD must be a string!');
end

% get radius
switch lower(method)
    case 'mean'
        r=(2*a+b)/3;
    case 'volumetric'
        r=(a*a*b)^(1/3);
    case 'authalic'
        a2=a^2; b2=b^2; sa2b2=sqrt(a2-b2);
        % avoid divide by zero for sphere
        if(sa2b2)
            r=sqrt((a2+(a*b2/sa2b2)*log((a+sa2b2)/b))/2);
        else
            r=a;
        end
    otherwise
        error('seizmo:mean_ellipsoid_radius:badMethod',...
            'Unknown METHOD: %s !',method);
end

end
