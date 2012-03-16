function [data]=makearrivals(data,phaselist,model)
%MAKEARRIVALS    Creates arrival info for SEIZMO data records
%
%    Usage:    data=makearrivals(data)
%              data=makearrivals(data,phaselist)
%              data=makearrivals(data,phaselist,model)
%
%    Description:
%     DATA=MAKEARRIVALS(DATA) adds seismic phase arrival time info to
%     records in SEIZMO struct DATA using the EV & ST header fields that
%     define the station/event geometry.  Uses the default model & phase
%     list from TAUPTIME.
%
%     DATA=MAKEARRIVALS(DATA,PHASELIST) uses the phase list defined by
%     PHASELIST instead of the TAUPTIME default.  
%
%     DATA=MAKEARRIVALS(DATA,PHASELIST,MODEL) uses the model defined by
%     MODEL instead of the TAUPTIME default.
%
%    Notes:
%
%    Examples:
%     % Use PREM and an expanded phase list:
%     data=makearrivals(data,'ttall','prem');
%
%     % Specify several S phases:
%     data=makearrivals(data,'sS,S,SS,SSS,SSSS,SKS,SKKS,ScS');
%
%    See also: FINDARRIVALS, ARRIVALS2PICKS, FINDPICKS

%     Version History:
%        Mar. 14, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 14, 2012 at 11:15 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% check data structure
data=checkheader(data,...
    'UNSET_ST_LATLON','ERROR',...
    'UNSET_EV_LATLON','ERROR',...
    'UNSET_ELEV','FIX',...
    'UNSET_DEPTH','FIX');

% default model/phaselist
if(nargin<2); phaselist=[]; end
if(nargin<3); model=[]; end

% get header info
[ev,st]=getheader(data,'ev','st');

% geographic to geocentric lat
ev(:,1)=geographic2geocentriclat(ev(:,1));
st(:,1)=geographic2geocentriclat(st(:,1));

% get arrivals for each record
for i=1:numel(data)
    data(i).misc.arrivals=tauptime('mod',model,'h',ev(i,4)/1000,...
        'ph',phaselist,'st',st(i,1:2),'ev',ev(i,1:2));
end

end
