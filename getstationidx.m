function [idx,stationcode]=getstationidx(data)
%GETSTATIONIDX    Returns index array for separating dataset into stations
%
%    Usage:    idx=getstationidx(data)
%              [idx,stationcode]=getstationidx(data)
%
%    Description: IDX=GETNETWOKIDX(DATA) returns an array of indices that
%     indicate records in SEIZMO structure DATA belonging to a station.  A
%     station is defined by the fields KNETWK and KSTNM.
%
%     [IDX,STREAMCODE]=GETSTATIONIDX(DATA) also returns the unique station
%     codes used to separate the stations.
%
%    Notes:
%     - Case insensitive; all characters are upper-cased.
%
%    Examples:
%     Break a dataset up into separate stations:
%      idx=getstationidx(data)
%      for i=1:max(idx)
%          stationdata{i}=data(idx==i);
%      end
%
%    See also: GETSTREAMIDX, GETNETWORKIDX, GETCOMPONENTIDX

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
[knetwk,kstnm]=getheader(data,'knetwk','kstnm');

% uppercase
knetwk=upper(knetwk);
kstnm=upper(kstnm);

% get station groups
[stationcode,idx,idx]=unique(strcat(knetwk,'.',kstnm));

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
