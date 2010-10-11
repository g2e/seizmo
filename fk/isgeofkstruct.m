function [lgc]=isgeofkstruct(fk)
%ISGEOFKSTRUCT    True if input is a fk struct as defined by GEOFK
%
%    Usage:    lgc=isgeofkstruct(geofk)
%
%    Description:
%     LGC=ISGEOFKSTRUCT(GEOFK) returns TRUE if GEOFK is a geofk struct or
%     FALSE otherwise.  GEOFK may be a nonscalar struct (LGC is still
%     either TRUE or FALSE).  See CHKGEOFKSTRUCT for more details.
%
%    Notes:
%
%    Examples:
%     % See if GEOFKXCVOLUME conforms:
%     isgeofkstruct(geofkxcvolume(xcdata,latlon,horzslow,freqrng))
%
%    See also: CHKFKSTRUCT, CHKFKARFSTRUCT, CHKGEOFKSTRUCT, ISFKARFSTRUCT,
%              CHKGEOFKARFSTRUCT, ISFKSTRUCT, ISGEOFKARFSTRUCT

%     Version History:
%        Oct. 10, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 10, 2010 at 14:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check struct
lgc=isempty(chkgeofkstruct(fk));

end
