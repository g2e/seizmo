function [data]=setevent(data,event)
%SETEVENT    Sets the event info for records in a SEIZMO dataset
%
%    Usage:    data=setevent(data,event)
%
%    Description: DATA=SETEVENT(DATA,EVENT) imports event info from a SOD
%     CSV struct (see READSODEVENTCSV) or a GLOBALCMT NDK struct (see
%     READNDK) into records in DATA.  EVENT must be a single event (ie the
%     struct must be a scalar).  Imported info is limited to time, location
%     and magnitude.
%
%    Notes:
%
%    Header changes: O, EVLA, EVLO, EVEL, EVDP, MAG, IMAGTYP, GCARC, AZ,
%     BAZ, DIST, KEVNM (for NDK events only)
%
%    Examples:
%     Import basic info from a quick CMT into some records:
%      ndk=readndk('quick.ndk');
%      data=setevent(data,ndk(33));
%
%    See also: READSODEVENTCSV, READNDK

%     Version History:
%        Dec.  1, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec.  1, 2009 at 22:45 GMT

% todo:

% check nargin
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data);
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% attempt header check
try
    % check header
    data=checkheader(data);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

% attempt rest
try
    % figure out if event is a sodcsv or ndk struct
    if(~isstruct(event) || ~isscalar(event))
        error('seizmo:setevent:badInput',...
            'EVENT must be a scalar struct!');
    end
    fields=fieldnames(event);
    sodfields={'time' 'latitude' 'longitude' ...
        'depth' 'magnitude' 'magnitudeType'};
    ndkfields={'year' 'month' 'day' 'hour' 'minute' 'seconds' ...
        'latitude' 'longitude' 'depth' 'mb' 'ms' 'name'};
    if(all(ismember(sodfields,fields)))
        % add info to header
        data=changeheader(data,'o 6utc',mat2cell(event.time,1),...
            'ev',[event.latitude event.longitude 0 event.depth*1000],...
            'mag',event.magnitude,'imagtyp',['i' event.magnitudeType]);
    elseif(all(ismember(ndkfields,fields)))
        % get the magnitude
        magtype='ims'; mag=event.ms;
        body=event.mb>event.ms;
        if(body); magtype='imb'; mag=event.mb; end
        % add info to header
        data=changeheader(data,...
            'o 6utc',{[event.year event.month event.day ...
            event.hour event.minute event.seconds]},...
            'ev',[event.latitude event.longitude 0 event.depth*1000],...
            'mag',mag,'imagtyp',magtype,'kevnm',event.name);
    else
        error('seizmo:setevent:badInput',...
            'EVENT struct type unknown!');
    end
    
    % update gcarc, az, baz, dist
    oldcheckheaderstate=get_checkheader_state;
    set_seizmocheck_state(true);
    data=checkheader(data,'all','ignore','old_delaz','fix');
    set_checkheader_state(oldcheckheaderstate);
    
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

end
