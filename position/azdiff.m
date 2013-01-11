function [d]=azdiff(az1,az2)
%AZDIFF    Returns the angle between azimuths
%
%    Usage:    d=azdiff(az1,az2)
%
%    Description:
%     D=AZDIFF(AZ1,AZ2) calculates the angle D from AZ1 to AZ2.  D is
%     always within +/-180.
%
%    Notes:
%
%    Examples:
%     % 340 --> 20 should give 40:
%     azdiff(340,20)
%
%     % 20 --> 340 should give -40:
%     azdiff(20,340)
%
%    See also: AZMEAN, SPHERICALINV, VINCENTYINV

%     Version History:
%        Feb.  9, 2012 - initial version
%        Oct. 10, 2012 - add AZMEAN to see also section
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 10, 2012 at 19:00 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% want angle d between 2 azimuths
d=mod(az2,360)-mod(az1,360);
d(d<-180)=d(d<-180)+360;
d(d>=180)=d(d>=180)-360;

end
