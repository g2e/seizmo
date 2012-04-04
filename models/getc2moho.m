function [moho]=getc2moho(lat,lon)
%GETC2MOHO    Returns Crust2.0 Moho depth at specified location(s)
%
%    Usage:    moho=getc2moho(lat,lon)
%
%    Description:
%     MOHO=GETC2MOHO(LAT,LON) returns the depth of the Moho discontinuity
%     at the positions specified.  The Moho is in kilometers depth from sea
%     level.  LAT & LON are expected to be scalars or equal sized arrays
%     and have units in degrees.  The depths are based on Crust2.0 by the
%     Scripps group.
%
%    Notes:
%
%    Examples:
%     % Plot a Crust2.0 moho depth map:
%     [x,y]=meshgrid(-179:2:179,89:-2:-89);
%     figure; imagesc(-179:2:179,89:-2:-89,getc2moho(y,x)); axis xy
%
%    See also: GETC2ELEV, GETC2THICK, GETCRUST2, CRUCOR, MANCOR

%     Version History:
%        May  17, 2010 - initial version
%        May  19, 2010 - name changed from GETMOHO
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 18:15 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check lat/lon
if(~isreal(lat) || ~isreal(lon))
    error('seizmo:getc2moho:badInput',...
        'LAT & LON must be real-valued arrays!');
elseif(~isequalsizeorscalar(lat,lon))
    error('seizmo:getc2moho:badInput',...
        'LAT & LON must be equal sized or scalar!');
end

% make the same size
[lat,lon]=expandscalars(lat,lon);

% make sure lat/lon values are in proper ranges
[lat,lon]=fixlatlon(lat,lon);

% load crust2.0
moho=load('crust2','moho');
moho=moho.moho;

% lat/lon to i/j to idx
i=fix((90-lat)/2)+1;
i(i==91)=90;
j=fix((lon+180)/2)+1;
j(j==181)=180;
idx=(j-1)*90+i;

% get moho for locations indicated
moho=moho(idx);

end
