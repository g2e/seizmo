function [paths]=getraypaths(phase,mod,evla,evlo,evdp,stla,stlo)
%GETRAYPATHS    Returns seismic phase paths for a set of stations
%
%    Usage:    paths=getraypaths(phase,mod,evla,evlo,evdp,stla,stlo)
%
%    Description:
%     PATHS=GETRAYPATHS(PHASE,MOD,EVLA,EVLO,EVDP,STLA,STLO) returns a
%     struct array containing phase paths for each PHASE/MODEL/EQ/STA pair.
%     PHASE should be a char or cellstr array and is case sensitive (uses
%     TauP to parse the phase name).  MOD must be a 1D model name
%     recognized by TauP.  Lat/Lon inputs must be in degrees.  Depth is in
%     kilometers (note that this is not like the SAC format)!
%
%     PATHS=GETRAYPATHS(DATA,PHASE,MOD) uses the station and event
%     positions in SEIZMO struct DATA to calculate the raypaths.
%
%    Notes:
%     - This just calls TAUPPATH in a loop.
%
%    Examples:
%     % Get some phase paths corresponding to a dataset:
%     paths=getraypaths(data,'P','prem');
%
%    See also: TAUPPATH, MANCOR, CRUST2LESS_RAYPATHS, TRIM_DEPTHS_RAYPATHS,
%              EXTRACT_UPSWING_RAYPATHS

%     Version History:
%        May  31, 2010 - initial version
%        June  3, 2010 - verbose support
%        June  4, 2010 - paths reshapes to match inputs
%        Mar.  7, 2011 - don't cover up error for debugging
%        May  19, 2011 - updated error message to reflect the usual issue
%        Feb. 27, 2012 - new usage form simplifies the typical case, return
%                        all paths
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 27, 2012 at 13:45 GMT

% todo:

% check nargin
error(nargchk(3,7,nargin));
if(~any(nargin==[3 7]))
    error('seizmo:getraypaths:badInput','Incorrect number of arguments!');
end

% handle data struct
if(isseizmo(phase))
    data=phase; phase=mod; mod=evla;
    data=checkheader(data,...
        'UNSET_ST_LATLON','ERROR',...
        'UNSET_EV_LATLON','ERROR',...
        'UNSET_ELEV','FIX',...
        'UNSET_DEPTH','FIX');
    [evla,evlo,evdp,stla,stlo]=getheader(data,...
        'evla','evlo','evdp','stla','stlo');
    evdp=evdp/1000;
end

% check inputs
if(ischar(phase)); phase=cellstr(phase); end
if(ischar(mod)); mod=cellstr(mod); end
if(~all(cellfun('isreal',{evla evlo evdp stla stlo})))
    error('seizmo:getraypaths:badLocation',...
        'All lat/lon/depth values must be real valued!');
elseif(~iscellstr(phase))
    error('seizmo:getraypaths:badPhase',...
        'PHASE must be char/cellstr array!');
elseif(~iscellstr(mod))
    error('seizmo:getraypaths:badModel',...
        'MOD must be char/cellstr array!');
elseif(~isequalsizeorscalar(phase,mod,evla,evlo,evdp,stla,stlo))
    error('seizmo:getraypaths:badSize',...
        'All inputs must be scalar or equal-sized!');
end

% expand scalars
[evla,evlo,evdp,stla,stlo,mod,phase]=expandscalars(...
    evla,evlo,evdp,stla,stlo,mod,phase);
nph=numel(evla);

% fix out-of-range lat/lon
[stla,stlo]=fixlatlon(stla,stlo);
[evla,evlo]=fixlatlon(evla,evlo);

% geographic to geocentric lat
evla=geographic2geocentriclat(evla);
stla=geographic2geocentriclat(stla);

% verbose
verbose=seizmoverbose;
if(verbose)
    disp('Getting Ray Path(s)');
    print_time_left(0,nph);
end

% loop over each set
try
    i=nph;
    tmp=tauppath('ph',phase{i},'mod',mod{i},'dep',evdp(i),...
        'ev',[evla(i) evlo(i)],'st',[stla(i) stlo(i)]);
    if(isempty(tmp))
        error('seizmo:getraypaths:badPath',...
            ['Could not retrieve path for EQ/STA pair: %d\n' ...
            'Maybe you do not have your ~/.taup file correct?'],i);
    end
    paths=tmp;
    if(verbose); print_time_left(1,nph); end
    for i=1:nph-1
        tmp=tauppath('ph',phase{i},'mod',mod{i},'dep',evdp(i),...
            'ev',[evla(i) evlo(i)],'st',[stla(i) stlo(i)]);
        if(isempty(tmp))
            error('seizmo:getraypaths:badPath',...
                'Could not retrieve path for EQ/STA pair: %d',i);
        end
        paths=[paths; tmp];
        if(verbose); print_time_left(i+1,nph); end
    end
catch
    error(lasterror);
end

end
