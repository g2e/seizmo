function [lgc]=isfkarfstruct(fk)
%ISFKARFSTRUCT    True if input is a fk struct as defined by FKARF
%
%    Usage:    lgc=isfkarfstruct(arf)
%
%    Description:
%     LGC=ISFKARFSTRUCT(ARF) returns TRUE if ARF is a fkarf struct or FALSE
%     otherwise.  ARF may be a nonscalar struct (LGC is still either TRUE
%     or FALSE).  See CHKFKARFSTRUCT for more details.
%
%    Notes:
%
%    Examples:
%     % See if FKARF produces a fkarf struct:
%     isfkarfstruct(fkarf(stla,stlo,50,201,0,0,[1/30 1/20]))
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
lgc=isempty(chkfkarfstruct(fk));

end
