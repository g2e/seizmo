function [lgc]=isstorms(storms)
%ISSTORMS    True if input is a storms struct
%
%    Usage:    lgc=isstorms(storms)
%
%    Description:
%     LGC=ISSTORMS(STORMS) returns TRUE if STORMS is a storms struct as
%     expected by MAPSTORMS and FALSE otherwise.
%
%    Notes:
%
%    Examples:
%     % See if the extratropical dataset conforms:
%     etstorms=load('etstorms.mat');
%     isstorms(etstorms)
%
%    See also: CHKSTORMS, MAPSTORMS, READ_HURDAT, READ_GISS_STORMDB

%     Version History:
%        Feb. 16, 2013 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 16, 2013 at 13:30 GMT

% todo:

error(nargchk(1,1,nargin));
lgc=isempty(chkstorms(storms));

end
