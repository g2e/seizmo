function [v]=smooth2d(v,nx,ny,edge)
%SMOOTH2D    2D data smoothing using a 2D Gaussian
%
%    Usage:    sv=smooth2d(v,n)
%              sv=smooth2d(v,nx,ny)
%              sv=smooth2d(v,nx,ny,edge)
%              sv=smooth2d(v,n,[],edge)
%
%    Description:
%     SV=SMOOTH2D(V,N) smooths the values in the 2D array V by convolving
%     it with a Gaussian with a characteristic falloff distance of N rows
%     or columns.  V is assumed to be a regularly spaced grid with equal
%     spacing in both the row and column directions.  N defines the
%     distance from the Gaussian's center (of unit amplitude) to the point
%     where the values have fallen off to 1/e.  Edge effects are accounted
%     for by using a equally smoothed unit surface for averaging.
%
%     SV=SMOOTH2D(V,NX,NY) allows for unequal spacing between the x
%     (column spacing) and the y (row spacing) directions.  For example,
%     the Gaussian may have a characteristic falloff distance of 4.5
%     x-spacings vs 3.2 y-spacings.
%
%     SV=SMOOTH2D(V,NX,NY,EDGE) or SV=SMOOTH2D(V,N,[],EDGE) indicates how
%     edge-effects are treated using the method EDGE.  EDGE must be a
%     string and can either be 'zeropad' or 'truncate'.  The zeropad option
%     just zero pads the array to be smoothed with zeros along the edges so
%     that the returned array has a tapered edge.  The truncate option will
%     truncate the gaussian smoother when it extends outside the array,
%     which makes the edges untapered but also under-smoothed.  The
%     default value of EDGE is 'truncate'.
%
%    Notes:
%
%    Examples:
%     % Compare PEAKS before/after applying a Gaussian with a falloff of 4
%     % grid spacings in both the x & y directions (0.5 units);
%     peaks;
%     surf(smooth2d(peaks,4));
%
%    See also: INTERP2, CONV2

%     Version History:
%        Feb.  4, 2010 - rewrite and added documentation
%        May  17, 2010 - updated documentation
%        July  9, 2010 - fixed nargchk
%        Apr.  2, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  2, 2012 at 17:10 GMT

% todo:

% check nargin
error(nargchk(2,4,nargin));

% check inputs
if(nargin<4); edge='truncate'; end
if(nargin==2 || isempty(ny)); ny=nx; end
if(ndims(v)~=2 || ~isnumeric(v))
    error('seizmo:smooth2d:badInput',...
        'V must be a 2D array of numeric values!');
elseif(~isscalar(nx) || ~isreal(nx) || nx<=0)
    error('seizmo:smooth2d:badInput',...
        'N or Nx must be a positive real scalar!');
elseif(~isscalar(ny) || ~isreal(ny) || ny<=0)
    error('seizmo:smooth2d:badInput',...
        'Ny must be a positive real scalar!');
elseif(~ischar(edge) || size(edge,1)~=1)
    error('seizmo:smooth2d:badInput',...
        'EDGE must be a string!');
end

% create gaussian
n=ceil(3*nx); gx=exp(-(abs(-n:n)/nx).^2);
n=ceil(3*ny); gy=exp(-(abs(-n:n)/ny).^2);

% smooth based on edge-effect method
switch lower(edge)
    case 'zeropad'
        % smooth array (tapers at edges due to zero-padding)
        gsum=sum(sum(gx'*gy));
        v=conv2(gy,gx,v,'same')/gsum;
    case 'truncate'
        % smooth array (accounts for zero-padding on edges)
        v=conv2(gy,gx,v,'same')./conv2(gy,gx,ones(size(v)),'same');
    otherwise
        error('seizmo:smooth2d:badInput',...
            'Unknown Edge-Effect Method: %s !',edge);
end

end
