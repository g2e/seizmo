function [paths,idx]=getraypaths(phaselist,model,evla,evlo,evdp,stla,stlo)
%GETRAYPATHS    Returns seismic phase paths for a set of stations
%
%    Usage:    paths=getraypaths(phaselist,model,evla,evlo,evdp,stla,stlo)
%              paths=getraypaths(data,phaselist,model)
%              [paths,idx]=getraypaths(...)
%
%    Description:
%     PATHS=GETRAYPATHS(PHASELIST,MODEL,EVLA,EVLO,EVDP,STLA,STLO) returns a
%     struct array containing phase paths for each PHASELIST/MODEL/EV/ST
%     set.  PHASELIST should be a char or cellstr array and is case
%     sensitive (uses TauP to parse the phase name).  MODEL must be a model
%     input recognized by TAUPPATH.  Lat/Lon inputs must be in degrees.
%     Depth is in kilometers (note that this is not the same as the output
%     from the EVDP header field)!
%
%     PATHS=GETRAYPATHS(DATA,PHASELIST,MODEL) uses the station and event
%     positions in SEIZMO struct DATA to calculate the raypaths.
%
%     [PATHS,IDX]=GETRAYPATHS(...) also returns the indices of the set to
%     which each path belongs.  This is mainly useful if using the DATA
%     input usage form.
%
%    Notes:
%     - This calls TAUPPATH in a loop.
%     - Latitudes are assumed to be geographic.
%
%    Examples:
%     % Get some phase paths corresponding to a dataset:
%     paths=getraypaths(data,'P','prem');
%
%    See also: TAUPPATH, MANCOR, CRUSTLESS_RAYPATHS, TRIM_DEPTHS_RAYPATHS,
%              EXTRACT_UPSWING_RAYPATHS, INSERT_DEPTHS_IN_RAYPATHS,
%              MAKEARRIVALS, FINDARRIVALS

%     Version History:
%        May  31, 2010 - initial version
%        June  3, 2010 - verbose support
%        June  4, 2010 - paths reshapes to match inputs
%        Mar.  7, 2011 - don't cover up error for debugging
%        May  19, 2011 - updated error message to reflect the usual issue
%        Feb. 27, 2012 - new usage form simplifies the typical case, return
%                        all paths
%        Mar. 14, 2012 - accept all mattaup model types, doc update
%        Aug.  6, 2012 - add latitude note to docs
%        Jan. 23, 2014 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 23, 2014 at 13:45 GMT

% todo:

% check nargin
error(nargchk(3,7,nargin));
if(~any(nargin==[3 7]))
    error('seizmo:getraypaths:badInput',...
        'Incorrect number of arguments!');
end

% handle data struct
if(isseizmo(phaselist))
    data=phaselist; phaselist=model; model=evla;
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
if(ischar(phaselist)); phaselist=cellstr(phaselist); end
if(ischar(model)); model=cellstr(model); end
if(~all(cellfun('isreal',{evla evlo evdp stla stlo})))
    error('seizmo:getraypaths:badLocation',...
        'All lat/lon/depth values must be real valued!');
elseif(~iscellstr(phaselist))
    error('seizmo:getraypaths:badPhase',...
        'PHASELIST must be char/cellstr array!');
elseif(~isequalsizeorscalar(phaselist,model,evla,evlo,evdp,stla,stlo))
    error('seizmo:getraypaths:badSize',...
        'All inputs must be scalar or equal-sized!');
end
if(~iscell(model)); model=num2cell(model); end

% expand scalars
[evla,evlo,evdp,stla,stlo,model,phaselist]=expandscalars(...
    evla,evlo,evdp,stla,stlo,model,phaselist);
nph=numel(evla);

% fix out-of-range lat/lon
[stla,stlo]=fixlatlon(stla,stlo);
[evla,evlo]=fixlatlon(evla,evlo);

% geographic to geocentric lat
evla=geographic2geocentriclat(evla);
stla=geographic2geocentriclat(stla);

% detail message
verbose=seizmoverbose;
if(verbose)
    disp('Getting Ray Path(s)');
    print_time_left(0,nph);
end

% now loop over each set
paths=[]; idx=[];
for i=1:nph
    tmp=tauppath('ph',phaselist{i},'mod',model{i},'dep',evdp(i),...
        'ev',[evla(i) evlo(i)],'st',[stla(i) stlo(i)]);
    if(isempty(tmp))
        error('seizmo:getraypaths:badPath',...
            'Could not retrieve path for EV/ST pair: %d',i);
    end
    if(isempty(paths)); paths=tmp; else paths=cat(1,paths,tmp); end
    idx=cat(1,idx,i*ones(numel(tmp),1));
    if(verbose); print_time_left(i,nph); end
end

end
