function [sod]=irisfetch2sod(event)
%IRISFETCH2SOD    Converts an irisFetch.Events struct to a SOD event struct
%
%    Usage:    sod=irisfetch2sod(ev)
%
%    Description:
%     SOD=IRISFETCH2SOD(EV) converts event info returned by IRISFETCH from
%     that struct format to one compatible with that used by SEIZMO's
%     READSODEVENTCSV, SOD2CMT, WRITESODEVENTCSV, and SETEVENT.
%
%    Notes:
%     - Only retains the preferred origin and magnitude info.  Additional
%       origins, magnitudes, and phase info is currently ignored.
%
%    Examples:
%     % Find GlobalCMTs that may be associated with events from an
%     % IRISFETCH search:
%     ev=irisFetch.Events('MinimumMagnitude',6.0,...
%         'boxcoordinates',[45 60 -150 -90]);
%     cmt=sod2cmt(irisfetch2sod(ev));
%
%     % One day's event list:
%     ev=irisFetch.Events('startTime','2011-01-07 00:00:00',...
%         'endTime','2011-01-08 00:00:00');
%     sod=irisfetch2sod(ev);
%     mmap('ev',[sod.latitude sod.longitude]);
%
%     % All events 10-20 degrees from the South Pole:
%     ev=irisFetch.Events('radialcoordinates',[-90 0 20 10]);
%     sod=irisfetch2sod(ev);
%     mmap('ev',[sod.latitude sod.longitude]);
%
%     % A large & deep event list:
%     ev=irisFetch.Events('minimumDepth',600,'minimumMagnitude',7);
%     sod=irisfetch2sod(ev);
%     figure; plot(sod.depth,sod.magnitude,'rp');
%
%     % Grab one specific event:
%     ev=irisFetch.Events('eventid','3028870');
%     sod=irisfetch2sod(ev)
%
%    See also: IRISFETCH, IRISFETCH2SEIZMO, SOD2CMT, WRITESODEVENTCSV,
%              READSODEVENTCSV, SETEVENT, SSIDX

%     Version History:
%        Feb.  9, 2014 - initial version
%        Feb. 20, 2014 - rename to irisfetch2sod
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 20, 2014 at 00:40 GMT

% todo:

% check number of inputs
error(nargchk(1,1,nargin));

% check input
if(~isstruct(event))
    error('seizmo:irisfetch2sod:badInput',...
        'EVENT must be a struct!');
elseif(isempty(event))
    error('seizmo:irisfetch2sod:badInput',...
        'EVENT is empty!');
elseif(any(~ismember({},fieldnames(event))))
    error('seizmo:irisfetch2sod:badInput',...
        'EVENT missing required struct fields!');
end

% number of events
nev=numel(event);

% translate to sod csv style struct
sod.time=datevec({event.PreferredTime}','yyyy-mm-dd HH:MM:SS.FFF');
sod.latitude=[event.PreferredLatitude]';
sod.longitude=[event.PreferredLongitude]';
sod.depth=[event.PreferredDepth]';
depthUnits={'kiloMETER'};
sod.depthUnits(1:nev,1)=depthUnits;
sod.magnitude=[event.PreferredMagnitudeValue]';
sod.magnitudeType={event.PreferredMagnitudeType}';
sod.catalog=getsubfield(event,'PreferredOrigin','Catalog')';
sod.contributor=getsubfield(event,'PreferredOrigin','Contributor')';
sod.type={event.Type}';
sod.name={event.PublicId}';
sod.flinnEngdahlRegionName={event.FlinnEngdahlRegionName}';
sod.flinnEngdahlRegionCode=[event.FlinnEngdahlRegionCode]';

end
