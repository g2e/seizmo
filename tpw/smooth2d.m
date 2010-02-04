function [v]=smooth2d(v,nx,ny)
%SMOOTH2D    2D data smoothing using a 2D Gaussian
%
%    Usage:    sv=smooth2d(v,n)
%              sv=smooth2d(v,nx,ny)
%
%    Description: SV=SMOOTH2D(V,N) smooths the values in the 2D array V
%     by convolving it with a Gaussian with a characteristic falloff of N.
%     V is assumed to be an evenly spaced grid with spacing equal in both
%     directions.  N is the distance (in grid spacings) from the Gaussian's
%     center to reach a value of 1/e.  Edge effects are accounted for by
%     using a equally smoothed unit surface for averaging.
%
%     SV=SMOOTH2D(V,NX,NY) allows for unequal spacing between the x
%     (column spacing) and the y (row spacing) directions.  For example,
%     the Gaussian may have a characteristic falloff distance of 4.5
%     x-spacings vs 3.2 y-spacings.
%
%    Notes:
%
%    Examples:
%     Compare PEAKS before/after applying a Gaussian with a falloff of 4
%     grid spacings in both the x & y directions (0.5 units);
%      peaks;
%      surf(smooth2d(peaks,4));
%
%    See also: INTERP2, CONV2

%     Version History:
%        Feb.  4, 2010 - rewrite and added documentation
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  4, 2010 at 17:10 GMT

% todo:

% check nargin
msg=nargchk(2,3,nargin);
if(~isempty(msg)); error(msg); end

% check inputs
if(nargin==2 || isempty(nx)); ny=nx; end
if(ndims(v)~=2 || ~isnumeric(v))
    error('seizmo:smooth2d:badInput',...
        'V must be a 2D array of numeric values!');
elseif(~isscalar(nx) || ~isreal(nx) || nx<=0)
    error('seizmo:smooth2d:badInput',...
        'N or Nx must be a positive real scalar!');
elseif(~isscalar(ny) || ~isreal(ny) || ny<=0)
    error('seizmo:smooth2d:badInput',...
        'Ny must be a positive real scalar!');
end

% create gaussian
n=ceil(3*nx); gx=exp(-(abs(-n:n)/nx).^2);
n=ceil(3*ny); gy=exp(-(abs(-n:n)/ny).^2);

% smooth array (accounts for zero-padding on edges)
v=conv2(gy,gx,v,'same')./conv2(gy,gx,ones(size(v)),'same');

end
