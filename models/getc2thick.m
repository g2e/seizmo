function [thick]=getc2thick(lat,lon)
%GETC2THICK    Returns Crust2.0 thickness at specified location(s)
%
%    Usage:    thick=getc2thick(lat,lon)
%
%    Description:
%     THICK=GETC2THICK(LAT,LON) returns the Crust2.0 crust thickness at the
%     positions specified.  Note that this has a 2x2 degree resolution.
%     LAT & LON are expected to be scalars or equal sized arrays and have
%     units in degrees.  Crust2.0 is by the Scripps group.
%
%    Notes:
%
%    Examples:
%     % Plot a Crust2.0 crustal thickness map:
%     [x,y]=meshgrid(-179:2:179,89:-2:-89);
%     figure; imagesc(-179:2:179,89:-2:-89,getc2thick(y,x)); axis xy
%
%    See also: GETC2MOHO, GETC2ELEV, GETCRUST2, CRUCOR, MANCOR

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
    error('seizmo:getc2thick:badInput',...
        'LAT & LON must be real-valued arrays!');
elseif(~isequalsizeorscalar(lat,lon))
    error('seizmo:getc2thick:badInput',...
        'LAT & LON must be equal sized or scalar!');
end

% make the same size
[lat,lon]=expandscalars(lat,lon);

% make sure lat/lon values are in proper ranges
[lat,lon]=fixlatlon(lat,lon);

% load crust2.0
thick=load('crust2','elev','moho');
thick=thick.elev/1000+thick.moho;

% lat/lon to i/j to idx
i=fix((90-lat)/2)+1;
i(i==91)=90;
j=fix((lon+180)/2)+1;
j(j==181)=180;
idx=(j-1)*90+i;

% get elev for locations indicated
thick=thick(idx);

end
