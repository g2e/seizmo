function [idx,streamcode]=getstreamidx(data)
%GETSTREAMIDX    Returns index array for separating dataset into streams
%
%    Usage:    idx=getstreamidx(data)
%              [idx,streamcode]=getstreamidx(data)
%
%    Description: IDX=GETSTREAMIDX(DATA) returns an array of indices that
%     indicate records in SEIZMO structure DATA belonging to a stream.  A
%     stream is defined by the fields KNETWK, KSTNM, KHOLE, and KCMPNM.
%     Only the third character of KCMPNM is allowed to vary between records
%     for a single stream.  This means that sites with multiple sensors
%     will be treated as separate streams as well as sensors that have
%     multiple streams of digitization (ie BH? and HH? channels would
%     be treated as separate streams).
%
%     [IDX,STREAMCODE]=GETSTREAMIDX(DATA) also returns the unique stream
%     codes used to separate the streams.
%
%    Notes:
%     - Case insensitive; all characters are upper-cased.
%
%    Examples:
%     Break a dataset up into separate streams:
%      idx=getstreamidx(data)
%      for i=1:max(idx)
%          streamdata{i}=data(idx==i);
%      end
%
%    See also: GETNETWORKIDX, GETSTATIONIDX, GETCOMPONENTIDX, GETSITEIDX

%     Version History:
%        June 28, 2009 - initial version
%        Jan. 29, 2010 - cleaned up unnecessary code
%        Aug. 11, 2010 - nargchk fix, doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 11, 2010 at 22:55 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% get header info
kname=upper(getheader(data,'kname'));

% truncate kcmpnm to first 2 characters
kname(:,4)=strnlen(kname(:,4),2);

% get stream groups
[streamcode,idx,idx]=unique(strcat(...
    kname(:,1),'.',kname(:,2),'.',kname(:,3),'.',kname(:,4)));

end
