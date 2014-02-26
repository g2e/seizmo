function [data]=irisfetch2seizmo(trace)
%IRISFETCH2SEIZMO    Converts an irisFetch.Traces struct to a SEIZMO struct
%
%    Usage:    data=irisfetch2seizmo(trace)
%
%    Description:
%     DATA=IRISFETCH2SEIZMO(TRACE) converts seismograms from IRISFETCH
%     from that struct format to the SEIZMO data struct format.  Values
%     are placed in their normal spots except a few fields that don't quite
%     fit:
%      TRACE.quality is set to the KDATRD field
%      TRACE.sensitivity is set to the SCALE field
%      TRACE.sensitivity* is placed in DATA.misc
%      TRACE.instrument is placed in DATA.misc
%      TRACE.sacpz stuff is in DATA.misc.sacpz (with some additional info)
%
%    Notes:
%     - IRISFETCH often returns an empty struct when there is actually data
%       at IRIS.  I have found that repeating the IRISFETCH command will
%       eventually return data so wrapping it in a while loop works.
%     - Timing from IRISFETCH is low resolution (millisecond accuracy) and
%       uses DATENUM to set TRACE.startTime & TRACE.endTime which can be
%       wrong when a leapsecond overlaps the time range of a record.
%
%    Examples:
%     % Download a few hours of vertical component data recorded by the
%     % Cameroon array (XB) during a particularly strong earthquake (this
%     % took about a minute to download on my home connection):
%     trace=[];
%     while(isempty(trace)) % force irisFetch.Traces to find something...
%         tic;
%         trace=irisFetch.Traces('XB','*','01','LHZ',...
%             '2006-01-31 19:15:51.590','2006-01-31 22:02:31.590',...
%             'includePZ');
%         toc;
%     end
%     data=irisfetch2seizmo(trace);
%     data=removetrend(synchronize(fix_cameroon(data),'b','first'));
%     plot0(squish(data,10),'normstyle','individual');
%
%     % Data from station CMB of the Berkeley network:
%     trace=[];
%     while(isempty(trace)) % force irisFetch.Traces to find something...
%         tic;
%         trace=irisFetch.Traces('BK','CMB','*','VHZ',...
%             '2006-01-31 19:15:51.590','2006-01-31 22:02:31.590',...
%             'includePZ');
%         toc;
%     end
%     data=irisfetch2seizmo(trace);
%     plot0(data);
%
%    See also: IRISFETCH, BSEIZMO, IRIS2SOD

%     Version History:
%        Feb.  7, 2014 - initial version
%        Feb.  8, 2014 - fixed cmpinc (0 is up, 90 horiz, 180 down)
%        Feb.  9, 2014 - doc update
%        Feb. 20, 2014 - rename to irisfetch2seizmo
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 20, 2014 at 00:40 GMT

% todo:

% check number of inputs
error(nargchk(1,1,nargin));

% check input
if(~isstruct(trace))
    error('seizmo:irisfetch2seizmo:badInput',...
        'TRACE must be a struct!');
elseif(isempty(trace))
    error('seizmo:irisfetch2seizmo:badInput',...
        'TRACE is empty!');
elseif(any(~ismember({'network' 'station' 'location' 'channel' ...
        'quality' 'latitude' 'longitude' 'elevation' 'depth' 'azimuth' ...
        'dip' 'data' 'sampleCount' 'sampleRate' 'startTime' 'sacpz' ...
        'sensitivity' 'sensitivityFrequency' 'instrument' ...
        'sensitivityUnits'},fieldnames(trace))))
    error('seizmo:irisfetch2seizmo:badInput',...
        'TRACE missing required struct fields!');
end

% number of records
nrecs=numel(trace);

% first preallocate output struct
verbose=seizmoverbose(false);
data(1:nrecs,1)=bseizmo([]);
seizmoverbose(verbose);

% populate data with records
[data.dep]=deal(trace.data);
[trace.data]=deal([]);

% dep*, sensitivity & sacpz
[depmin,depmen,depmax]=deal(nan(nrecs,1));
b=serial2gregorian([trace.startTime]','doytime');
e=serial2gregorian([trace.endTime]','doytime');
name=strcat('SAC_PZs_',{trace.network}','_',{trace.station}','_',...
    {trace.channel}','_',{trace.location}','_',...
    datestr([trace.startTime]','yyyy.mm.dd.HH.MM.SS.FFF'),'_',...
    datestr([trace.endTime]','yyyy.mm.dd.HH.MM.SS.FFF'));
for i=1:nrecs
    depmen(i)=nanmean(data(i).dep(:));
    depmin(i)=min(data(i).dep(:));
    depmax(i)=max(data(i).dep(:));
    data(i).misc.instrument=trace(i).instrument;
    data(i).misc.sensitivity=trace(i).sensitivity;
    data(i).misc.sensitivityFrequency=trace(i).sensitivityFrequency;
    data(i).misc.sensitivityUnits=trace(i).sensitivityUnits;
    if(isempty(trace(i).sacpz.constant))
        data(i).misc.has_sacpz=false;
        data(i).misc.sacpz=[];
    else
        data(i).misc.has_sacpz=true;
        data(i).misc.sacpz.path='./';
        data(i).misc.sacpz.name=name{i};
        data(i).misc.sacpz.knetwk=trace(i).network;
        data(i).misc.sacpz.kstnm=trace(i).station;
        data(i).misc.sacpz.khole=trace(i).location;
        data(i).misc.sacpz.kcmpnm=trace(i).channel;
        data(i).misc.sacpz.b=b(i,:); % we know it is valid from this time
        data(i).misc.sacpz.e=e(i,:); % to this time ...
        data(i).misc.sacpz.z=trace(i).sacpz.zeros;
        data(i).misc.sacpz.p=trace(i).sacpz.poles;
        data(i).misc.sacpz.k=trace(i).sacpz.constant;
    end
end

% insert info into the header
data=changeheader(data,'npts',[trace.sampleCount],...
    'delta',1./[trace.sampleRate],'z',b,'b',0,'iztype','ib',...
    'cmpaz',[trace.azimuth],'cmpinc',90+[trace.dip],...
    'stla',[trace.latitude],'stlo',[trace.longitude],...
    'stel',[trace.elevation],'stdp',[trace.depth],...
    'knetwk',{trace.network},'kstnm',{trace.station},...
    'khole',{trace.location},'kcmpnm',{trace.channel},...
    'kdatrd',{trace.quality},'scale',[trace.sensitivity]);

% give them filenames
data=genname(data);

end
