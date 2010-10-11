function [lgc]=isgeofkarfstruct(fk)
%ISGEOFKARFSTRUCT    True if input is a fk struct as defined by GEOFKARF
%
%    Usage:    lgc=isgeofkarfstruct(arf)
%
%    Description:
%     LGC=ISGEOFKARFSTRUCT(ARF) returns TRUE if ARF is a geofkarf struct or
%     FALSE otherwise.  ARF may be a nonscalar struct (LGC is still either
%     TRUE or FALSE).  See CHKGEOFKARFSTRUCT for more details.
%
%    Notes:
%
%    Examples:
%     % See if GEOFKARF produces a geofkarf struct:
%     [lat,lon]=meshgrid(-20:20,-20:20);
%     isgeofkarfstruct(geofkarf([stla stlo],[lat(:) lon(:)],...
%         27.5:.5:32.5,[5 10],30,1/26.3))
%
%    See also: CHKFKSTRUCT, CHKFKARFSTRUCT, CHKGEOFKSTRUCT, ISFKSTRUCT,
%              CHKGEOFKARFSTRUCT, ISGEOFKSTRUCT, ISGEOFKARFSTRUCT

%     Version History:
%        Oct. 10, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 10, 2010 at 14:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check struct
lgc=isempty(chkgeofkarfstruct(fk));

end
