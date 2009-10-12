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
[knetwk,kstnm,khole,kcmpnm]=...
    getheader(data,'knetwk','kstnm','khole','kcmpnm');

% uppercase
knetwk=upper(knetwk);
kstnm=upper(kstnm);
khole=upper(khole);
kcmpnm=upper(kcmpnm);

% get station groups
[componentcode,idx,idx]=...
    unique(strcat(knetwk,'.',kstnm,'.',khole,'.',kcmpnm));

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end
