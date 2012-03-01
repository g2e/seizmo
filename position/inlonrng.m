function [tf]=inlonrng(rng,lon)
%INLONRNG    Returns TRUE for longitudes within the specified range
%
%    Usage:    tf=inlonrng(rng,lon)
%
%    Description:
%     TF=INLONRNG(RNG,LON) returns TRUE or FALSE depending on if LON is in
%     the longitude range specified by RNG.  RNG must be [MINLON MAXLON]
%     and handles wraparound of longitudes.
%
%    Notes:
%
%    Examples:
%     % A few tests that should return TRUE:
%     inlonrng([-10 10],360)
%     inlonrng([170 190],-180)
%     inlonrng([170 190],180)
%     inlonrng([-190 -170],-180)
%     inlonrng([-190 -170],180)
%     inlonrng([350 370],0)
%
%     % A few tougher cases that should return TRUE:
%     inlonrng([0 0],360)
%     inlonrng([180 180],-180)
%     inlonrng([180 180],180)
%     inlonrng([-180 -180],-180)
%     inlonrng([-180 -180],180)
%     inlonrng([360 360],0)
%
%     % Yet more cases that should return TRUE:
%     inlonrng([0 360],0)
%     inlonrng([0 360],360)
%     inlonrng([-180 180],-180)
%     inlonrng([-180 180],180)
%
%    See also: LONMOD, LATMOD, FIXLATLON

%     Version History:
%        Feb. 29, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 29, 2012 at 17:25 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check args
if(~isnumeric(rng) || ~isreal(rng) || size(rng,2)~=2 || ndims(rng)>2)
    error('seizmo:inlonrng:badInput',...
        'RNG must be a real-valued array given as [MINLON MAXLON]!');
elseif(any(rng(:,1)>rng(:,2)))
    error('seizmo:inlonrng:badInput',...
        'RNG must be a real-valued array given as [MINLON MAXLON]!');
elseif(~isnumeric(lon) || ~isreal(lon))
    error('seizmo:inlonrng:badInput',...
        'LON must be a real-valued array!');
elseif(size(rng,1)~=1 && ~isscalar(lon) && numel(lon)~=size(rng,1))
    error('seizmo:inlonrng:badInput',...
        'Number of rows in RNG must match the number of elements in LON!');
end

% slow but infallible
% - there might be a better way with azdiff...
drng=diff(rng,1,2);
rng(:,1)=lonmod(rng(:,1));
rng(:,2)=rng(:,1)+drng;
lon=lonmod(lon);
tf=(lon>=rng(:,1) & lon<=rng(:,2)) ...
    | (lon+360>=rng(:,1) & lon+360<=rng(:,2)) ...
    | (lon-360>=rng(:,1) & lon-360<=rng(:,2));

end
