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
%                    BAZ, DIST, KEVNM
%
%    Examples:
%     Import basic info from a quick CMT into some records:
%      cmt=load('globalcmt_quick');
%      data=setevent(data,ssidx(cmt,33));
%
%    See also: READSODEVENTCSV, READNDK, SSIDX

%     Version History:
%        Dec.  1, 2009 - initial version
%        Jan. 26, 2010 - seizmoverbose support
%        Jan. 30, 2010 - fixed bug in call to CHECKHEADER
%        Mar. 30, 2010 - imagtyp for sod is optional, kevnm/name for sod is
%                        new optional, use mw for ndk
%        July 30, 2010 - minor touches for new ndk setup
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 30, 2010 at 22:25 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check data structure & header
data=checkheader(data);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt rest
try
    % verbosity
    verbose=seizmoverbose;

    % number of records
    nrecs=numel(data);
    
    % figure out if event is a sodcsv or ndk struct
    if(~isstruct(event) || ~isscalar(event))
        error('seizmo:setevent:badInput',...
            'EVENT must be a scalar struct!');
    end
    fields=fieldnames(event);
    sodfields={'time' 'latitude' 'longitude' 'depth' 'magnitude'};
    ndkfields={'year' 'month' 'day' 'hour' 'minute' 'seconds' ...
        'latitude' 'longitude' 'depth' 'scalarmoment' 'exponent' 'name'};
    if(all(ismember(sodfields,fields)))
        % detail message
        if(verbose)
            disp('Importing Event Info');
            print_time_left(0,nrecs);
        end
        
        % check for magtype
        imagtyp=nan;
        if(ismember('magnitudeType',fields))
            imagtyp=['i' event.magnitudeType];
        end
        
        % check for name
        kevnm=nan;
        if(ismember('name',fields))
            kevnm=event.name;
        end
        
        % add info to header
        data=changeheader(data,'o 6utc',mat2cell(event.time,1),...
            'ev',[event.latitude event.longitude 0 event.depth*1000],...
            'mag',event.magnitude,'imagtyp',imagtyp,'kevnm',kevnm);
    elseif(all(ismember(ndkfields,fields)))
        % detail message
        if(verbose)
            disp('Importing Event Info');
            print_time_left(0,nrecs);
        end
        
        % get the magnitude
        mw=(2/3)*log10(event.scalarmoment*10^event.exponent)-10.7;
        mw=round(mw*10)/10;
        
        % add info to header
        data=changeheader(data,...
            'o 6utc',{[event.year event.month event.day ...
            event.hour event.minute event.seconds]},...
            'ev',[event.latitude event.longitude 0 event.depth*1000],...
            'mag',mw,'imagtyp','imw','kevnm',event.name);
    else
        error('seizmo:setevent:badInput',...
            'EVENT struct type unknown!');
    end
    
    % update gcarc, az, baz, dist
    oldcheckheaderstate=checkheader_state;
    checkheader_state(true);
    data=checkheader(data,'all','ignore','old_delaz','fix');
    checkheader_state(oldcheckheaderstate);
    
    % detail message
    if(verbose); print_time_left(nrecs,nrecs); end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
