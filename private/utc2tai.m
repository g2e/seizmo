function [tai]=utc2tai(utc)
%UTC2TAI    Convert UTC time to TAI time
%
%    Description: UTC2TAI(UTC) returns
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Usage:    tai=utc2tai(utc)
%
%    Examples:
%
%    See also: tai2utc, cleanutc

%     Version History:
%        Nov.  1, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  1, 2008 at 17:40 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% require Nx5 or Nx6 numeric
sz=size(utc);
ndates=prod(sz)/sz(2);
if(~isnumeric(utc) || ~any(sz(2)==[5 6]))
    error('SAClab:utc2tai:badInput','UTC array must be Nx5 or Nx6!');
end

% permute for easier expansion and then make 2D
utc=permute(utc,[2 1 3:numel(sz)]);
utc=utc(:,:).';

% fix minutes and hours
utc(:,end-2)=utc(:,end-2)+floor(utc(:,end-1)/60);
utc(:,end-1)=mod(utc(:,end-1),60);
utc(:,end-3)=utc(:,end-3)+floor(utc(:,end-2)/24);
utc(:,end-2)=mod(utc(:,end-2),24);

% clean up date (datenum cannot handle low months)
if(sz(2)==5); utc=[doy2cal(fixdates(utc(:,1:2))) utc(:,3:5)];
else utc(:,1:3)=fixdates(utc(:,1:3));
end

% TAI time of start of minute
tai=datevec(datenum([utc(:,1:5) zeros(ndates,1)]));
tai(:,6)=tai(:,6)+totalleaps(utc(:,1:3))+10;

% TAI time
tai(:,6)=tai(:,6)+utc(:,6);

% fix TAI time
tai(:,5)=tai(:,5)+floor(tai(:,6)/60);
tai(:,6)=mod(tai(:,6),60);
tai(:,4)=tai(:,4)+floor(tai(:,5)/60);
tai(:,5)=mod(tai(:,5),60);
tai(:,3)=tai(:,3)+floor(tai(:,4)/24);
tai(:,4)=mod(tai(:,4),24);
tai(:,1:3)=fixdates(tai(:,1:3));

% reshape and convert to doy if necessary
if(sz(2)==5); tai=[cal2doy(tai(:,1:3)) tai(:,4:6)]; end
tai=reshape(tai,sz);

end
