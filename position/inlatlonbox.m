function [tf]=inlatlonbox(latrng,lonrng,latlon)
%INLATLONBOX    Returns TRUE for positions within the specified region
%
%    Usage:    tf=inlatlonbox(latrng,lonrng,latlon)
%
%    Description:
%     TF=INLATLONBOX(LATRNG,LONRNG,LATLON) returns TRUE or FALSE depending
%     on if the positions given by LATLON are within the region specified
%     by LATRNG & LONRNG.  LATRNG & LONRNG should be specified as [MIN MAX]
%     and LATLON should be given as [LAT LON].  TF is a NROWSx1 logical
%     array where NROWS is the number of rows in LATLON.
%
%    Notes:
%
%    Examples:
%     % Make a region covering 1/4th the Earth and see how close to 1/4th
%     % of a random set of positions is in the region:
%     sum(inlatlonbox([0 90],[0 180],randlatlon(1000)))/1000
%
%    See also: INLONRNG, FIXLATLON

%     Version History:
%        Feb. 29, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 29, 2012 at 17:25 GMT

% todo:

% check nargin
error(nargchk(3,3,nargin));

% check args
if(~isnumeric(latrng) || ~isreal(latrng) ...
        || size(latrng,2)~=2 || ndims(latrng)>2 ...
        || any(latrng(:,1)>latrng(:,2)))
    error('seizmo:inlonrng:badInput',...
        'LATRNG must be a real-valued array given as [MINLAT MAXLAT]!');
elseif(~isnumeric(lonrng) || ~isreal(lonrng) ...
        || size(lonrng,2)~=2 || ndims(lonrng)>2 ...
        || any(lonrng(:,1)>lonrng(:,2)))
    error('seizmo:inlonrng:badInput',...
        'LONRNG must be a real-valued array given as [MINLON MAXLON]!');
elseif(~isnumeric(latlon) || ~isreal(latlon) ...
        || size(latlon,2)~=2 || ndims(latlon)>2)
    error('seizmo:inlonrng:badInput',...
        'LATLON must be a real-valued array given as [LAT LON]!');
end
sz=[size(latrng,1) size(lonrng,1) size(latlon,1)];
if(sum(sz~=1)>1 && numel(unique(sz(sz~=1)))>1)
    error('seizmo:inlonrng:badInput',...
        'All inputs must have the 1 row or the same number of rows!');
end

% find those in box
tf=latlon(:,1)>=latrng(:,1) & latlon(:,1)<=latrng(:,2) ...
    & inlonrng(lonrng,latlon(:,2));

end
