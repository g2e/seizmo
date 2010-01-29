function [idx,componentcode]=getcomponentidx(data)
%GETCOMPONENTIDX    Returns index array to separate dataset into components
%
%    Usage:    idx=getcomponentidx(data)
%              [idx,componentcode]=getcomponentidx(data)
%
%    Description: IDX=GETCOMPONENTIDX(DATA) returns an index array that
%     indicates records in SEIZMO structure DATA belonging to a component.
%     A component is defined by the fields KNETWK, KSTNM, KHOLE, & KCMPNM.
%     This means that sites with multiple sensors will be separated as well
%     as multiple streams of digitization (ie BH? and HH? streams would
%     be treated as separate components).
%
%     [IDX,COMPONENTCODE]=GETCOMPONENTIDX(DATA) also returns the unique
%     component codes used to separate the components.
%
%    Notes:
%     - Case insensitive; all characters are upper-cased.
%
%    Examples:
%     Break a dataset up into separate components:
%      idx=getcomponentidx(data)
%      for i=1:max(idx)
%          componentdata{i}=data(idx==i);
%      end
%
%    See also: GETNETWORKIDX, GETSTATIONIDX, GETSTREAMIDX

%     Version History:
%        June 28, 2009 - initial version
%        Jan. 29, 2010 - cleaned up unnecessary code
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 29, 2010 at 23:20 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% get header info
kname=upper(getheader(data,'kname'));

% get component groups
[componentcode,idx,idx]=unique(strcat(...
    kname(:,1),'.',kname(:,2),'.',kname(:,3),'.',kname(:,4)));

end
