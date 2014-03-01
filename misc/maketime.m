function [t]=maketime(b,npts,delta)
%MAKETIME    Makes discrete time values given B, NPTS & DELTA
%
%    Usage:    t=maketime(b,npts,delta)
%
%    Description:
%     T=MAKETIME(B,NPTS,DELTA) returns the discrete time points for a given
%     starting sample time B, number of points NPTS, and sample spacing
%     DELTA.  B & DELTA must be in the same units (e.g., seconds).  All
%     inputs must be scalar.  T is a NPTSx1 vector.
%
%    Notes:
%     - MAKETIME uses the colon operator to limit floating point issues.
%       Note that the following two formulas for T will be different at
%       floating point precision:
%        T=B+(0:NPTS-1)*DELTA
%        T=B+(0:DELTA:DELTA*(NPTS-1))
%       MAKETIME uses the later.
%
%    Examples:
%     % The following is basically the equivalent
%     % of SEIZMO's MAKEUNEVEN function on 1 record:
%     [b,npts,delta]=getheader(data(1),'b','npts','delta');
%     data(1).ind=maketime(b,npts,delta);
%     data(1)=changeheader(data(1),'leven',false);
%
%     % LINSPACE, the colon (:) operator, & MAKETIME are quite similar:
%     b=0; npts=5; delta=1; e=b+(npts-1)*delta;
%     linspace(b,e,npts)
%     maketime(b,npts,delta)
%     b:delta:e
%
%    See also: MAKEUNEVEN

%     Version History:
%        Feb. 20, 2014 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 20, 2014 at 11:15 GMT

% todo:

% check number of inputs
error(nargchk(3,3,nargin));

% check inputs
if(~isscalar(b) || ~isnumeric(b))
    error('seizmo:maketime:badInput',...
        'B must be a scalar number!');
elseif(~isscalar(npts) || ~isreal(npts) || npts<0 || npts~=fix(npts))
    error('seizmo:maketime:badInput',...
        'NPTS must be a positive integer!');
elseif(~isscalar(delta) || ~isnumeric(delta))
    error('seizmo:maketime:badInput',...
        'DELTA must be a scalar number!');
end

% get discrete times
% - note that the delta is in the : operation
%   as this does better for floating point
t=b+(0:delta:delta*(npts-1));

end
