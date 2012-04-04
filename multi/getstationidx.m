function [idx,stationcode]=getstationidx(data)
%GETSTATIONIDX    Returns index array for separating dataset into stations
%
%    Usage:    idx=getstationidx(data)
%              [idx,stationcode]=getstationidx(data)
%
%    Description:
%     IDX=GETNETWOKIDX(DATA) returns an array of indices that indicate
%     records in SEIZMO structure DATA belonging to a station.  A station
%     is defined by the fields KNETWK and KSTNM.
%
%     [IDX,STREAMCODE]=GETSTATIONIDX(DATA) also returns the unique station
%     codes used to separate the stations.
%
%    Notes:
%     - Case insensitive; all characters are upper-cased.
%
%    Examples:
%     % Break a dataset up into separate stations:
%     idx=getstationidx(data)
%     for i=1:max(idx)
%         stationdata{i}=data(idx==i);
%     end
%
%    See also: GETSTREAMIDX, GETNETWORKIDX, GETCOMPONENTIDX, GETSITEIDX

%     Version History:
%        June 28, 2009 - initial version
%        Jan. 29, 2010 - cleaned up unnecessary code
%        Aug. 11, 2010 - nargchk fix, doc update
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 23:15 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% get header info
[knetwk,kstnm]=getheader(data,'knetwk','kstnm');

% uppercase
knetwk=upper(knetwk);
kstnm=upper(kstnm);

% get station groups
[stationcode,idx,idx]=unique(strcat(knetwk,'.',kstnm));

end
