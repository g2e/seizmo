function [idx,sitecode]=getsiteidx(data)
%GETSITEIDX    Returns index array for separating dataset by sites
%
%    Usage:    idx=getsiteidx(data)
%              [idx,streamcode]=getsiteidx(data)
%
%    Description:
%     IDX=GETSITEIDX(DATA) returns an array of indices that indicate
%     records in SEIZMO structure DATA belonging to a specific site (aka
%     location).  A site is defined by the fields KNETWK, KSTNM, and KHOLE.
%     Usually this will separate records into groups based on a specific
%     sensor or sample rate.  The conventions employed vary from network to
%     network unfortunately so this isn't perfect by any means.
%
%     [IDX,STREAMCODE]=GETSITEIDX(DATA) also returns the unique site
%     codes used to separate the sites.
%
%    Notes:
%     - Case insensitive; all characters are upper-cased.
%
%    Examples:
%     % Break a dataset up into separate sites:
%     idx=getsiteidx(data)
%     for i=1:max(idx)
%         sitedata{i}=data(idx==i);
%     end
%
%    See also: GETNETWORKIDX, GETSTATIONIDX, GETCOMPONENTIDX, GETSTREAMIDX

%     Version History:
%        Aug. 11, 2010 - initial version
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 22:55 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% get header info
kname=upper(getheader(data,'kname'));

% get stream groups
[sitecode,idx,idx]=unique(strcat(...
    kname(:,1),'.',kname(:,2),'.',kname(:,3)));

end
