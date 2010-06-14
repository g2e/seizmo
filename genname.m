function [data]=genname(data)
%GENNAME    Generic filenames for SEIZMO records
%
%    Usage:    data=genname(data)
%
%    Description: DATA=GENNAME(DATA) generates a generic filename for each
%     record in SEIZMO struct DATA using the KNAME header field group.
%     Care is taken to avoid duplicate filenames in the same dataset but
%     not with existing files or other datasets.  This is useful for simple
%     and short filenames that indicate the station the record corresponds
%     to.  The .name field of each record is replaced by:
%
%         SAC.KNETWK.KSTNM.KHOLE.KCMPNM.IDX
%
%     where IDX is a 2+ digit number starting at 00 and incremented for
%     records with the same KNAME header field set.
%
%    Notes:
%     - this matches with my perl code naming schemes
%
%    Examples:
%     Read, set the name, write:
%      writeseizmo(genname(readseizmo('*')));
%
%    See also: CHANGENAME

%     Version History:
%        June 13, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 13, 2010 at 01:45 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% get necessary info
kname=getheader(data,'kname');

% combine into a simple string
kname=upper(kname);
kname=strcat(kname(:,1),'.',kname(:,2),'.',kname(:,3),'.',kname(:,4));

% pass through unique for finding repeats
[uname,idx,idx]=unique(kname);

% occurance index (starts with zero)
oi=nan(size(idx));
for i=1:numel(uname)
    cmps=idx==i;
    oi(cmps)=0:sum(cmps)-1;
end

% now finish off the string
digits=num2str(max(2,log10(max(oi))+1)); % at least 2 digits
kname=strcat('SAC.',kname,'.',num2str(oi,['%0' digits 'd']));

% rename
[data.name]=deal(kname{:});

end
