function [lgc]=isfkstruct(fk)
%ISFKSTRUCT    True if input is a fk struct as defined by FKMAP/VOLUME/4D
%
%    Usage:    lgc=isfkstruct(fk)
%
%    Description:
%     LGC=ISFKSTRUCT(FK) returns TRUE if FK is a fk struct or FALSE
%     otherwise.  FK may be a nonscalar struct (LGC is still either TRUE or
%     FALSE).  See CHKFKSTRUCT for more details.
%
%    Notes:
%
%    Examples:
%     % See if FKMAP produces a fk struct:
%     isfkstruct(fkmap(data,50,201,[.04 .05]))
%
%    See also: CHKFKSTRUCT, CHKFKARFSTRUCT, CHKGEOFKSTRUCT, ISFKARFSTRUCT,
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
lgc=isempty(chkfkstruct(fk));

end
