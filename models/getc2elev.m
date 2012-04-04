function [elev]=getc2elev(lat,lon)
%GETC2ELEV    Returns Crust2.0 elevation at specified location(s)
%
%    Usage:    elev=getc2elev(lat,lon)
%
%    Description:
%     ELEV=GETC2ELEV(LAT,LON) returns Crust2.0 elevations at the positions
%     specified.  Note that this has a 2x2 degree resolution.  LAT & LON
%     are expected to be scalars or equal sized arrays and have units in
%     degrees.  Crust2.0 is by the Scripps group.
%
%    Notes:
%
%    Examples:
%     % Plot a Crust2.0 elevation map:
%     [x,y]=meshgrid(-179:2:179,89:-2:-89);
%     figure; imagesc(-179:2:179,89:-2:-89,getc2elev(y,x)); axis xy
%
%    See also: GETC2MOHO, GETC2THICK, GETCRUST2, CRUCOR, MANCOR

%     Version History:
%        May  19, 2010 - initial version
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 17:55 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check lat/lon
if(~isreal(lat) || ~isreal(lon))
    error('seizmo:getc2elev:badInput',...
        'LAT & LON must be real-valued arrays!');
elseif(~isequalsizeorscalar(lat,lon))
    error('seizmo:getc2elev:badInput',...
        'LAT & LON must be equal sized or scalar!');
end

% make the same size
[lat,lon]=expandscalars(lat,lon);

% make sure lat/lon values are in proper ranges
[lat,lon]=fixlatlon(lat,lon);

% load crust2.0
elev=load('crust2','elev');
elev=elev.elev;

% lat/lon to i/j to idx
i=fix((90-lat)/2)+1;
i(i==91)=90;
j=fix((lon+180)/2)+1;
j(j==181)=180;
idx=(j-1)*90+i;

% get elev for locations indicated
elev=elev(idx);

end
