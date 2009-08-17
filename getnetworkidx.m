function [idx,networkcode]=getnetworkidx(data)
%GETNETWORKIDX    Returns index array for separating dataset into networks
%
%    Usage:    idx=getnetworkidx(data)
%              [idx,networkcode]=getnetworkidx(data)
%
%    Description: IDX=GETNETWOKIDX(DATA) returns an array of indices that
%     indicate records in SEIZMO structure DATA belonging to a network.  A
%     network is defined by the field KNETWK.
%
%     [IDX,STREAMCODE]=GETNETWORKIDX(DATA) also returns the unique network
%     codes used to separate the networks.
%
%    Notes:
%     - Case insensitive; all characters are upper-cased.
%
%    Examples:
%     Break a dataset up into separate networks:
%      idx=getnetworkidx(data)
%      for i=1:max(idx)
%          networkdata{i}=data(idx==i);
%      end
%
%    See also: getstreamidx, getstationidx, getcomponentidx

%     Version History:
%        June 28, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:05 GMT

% todo:

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data);
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% get header info
[knetwk]=getheader(data,'knetwk');

% uppercase
knetwk=upper(knetwk);

% get station groups
[networkcode,idx,idx]=unique(knetwk);

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
