function [varargout]=mapstations(data,varargin)
%MAPSTATIONS    Map station/earthquake locations of SEIZMO records
%
%    Usage:    mapstations(data)
%              mapstations(...,'stations',[lat lon],...)
%              mapstations(...,'stationmarker',symstr,...)
%              mapstations(...,'stationmarkersize',val,...)
%              mapstations(...,'events',[lat lon],...)
%              mapstations(...,'eventmarker',symstr,...)
%              mapstations(...,'eventmarkersize',symstr,...)
%              mapstations(...,'gshhs',res,...)
%              mapstations(...,'proj',proj,...)
%              mapstations(...,'projopt',{'opt',val,...},...)
%              mapstations(...,'gridopt',{'opt',val,...},...)
%              mapstations(...,'fgcolor',color,...)
%              mapstations(...,'bgcolor',color,...)
%              mapstations(...,'sea',color,...)
%              mapstations(...,'land',color,...)
%              mapstations(...,'border',color,...)
%              mapstations(...,'axis',ax,...)
%              ax=mapstations(...)
%
%    Description: MAPSTATIONS(DATA) creates a map showing the station and
%     earthquake locations stored in the headers of records of SEIZMO
%     struct DATA.  The map is a global map using the Robinson projection.
%     Stations are plotted as yellow circles and events are plotted as
%     5-pointed stars.
%
%     MAPSTATIONS(...,'STATIONS',[LAT LON],...)
%     MAPSTATIONS(...,'STATIONMARKER',SYMSTR,...)
%     MAPSTATIONS(...,'STATIONMARKERSIZE',VAL,...)
%     MAPSTATIONS(...,'EVENTS',[LAT LON],...)
%     MAPSTATIONS(...,'EVENTMARKER',SYMSTR,...)
%     MAPSTATIONS(...,'EVENTMARKERSIZE',SYMSTR,...)
%     MAPSTATIONS(...,'GSHHS',RES,...)
%     MAPSTATIONS(...,'PROJ',PROJ,...)
%     MAPSTATIONS(...,'PROJOPT',{'OPT',VAL,...},...)
%     MAPSTATIONS(...,'GRIDOPT',{'OPT',VAL,...},...)
%     MAPSTATIONS(...,'FGCOLOR',COLOR,...)
%     MAPSTATIONS(...,'BGCOLOR',COLOR,...)
%     MAPSTATIONS(...,'SEA',COLOR,...)
%     MAPSTATIONS(...,'LAND',COLOR,...)
%     MAPSTATIONS(...,'BORDER',COLOR,...)
%     MAPSTATIONS(...,'AXIS',AX,...)
%     AX=MAPSTATIONS(...)
%      see MAPLOCATIONS for details on these options.  Note that using the
%      STATIONS or EVENTS option overrides the info found in DATA.
%
%    Notes:
%
%    Examples:
%     Show locations of stations in a dataset:
%      mapstations(data);
%
%     This replaces the old MAPSTATIONS2:
%      ax=mapstations(data);
%      ev=getheader(data(1),'ev');
%      mapeventgrid(ax,ev(1),ev(2));
%
%     This replaces the old MAPCLUSTERS:
%      ax=mapstations(data);
%      h=findobj(ax,'tag','stations');
%      set(h,'cdata',grp.color(grp.T,:));
%
%     Show Africa stations in a map with fancy border:
%      mapstations(data,'po',{'lat',[-40 40],'lon',[-30 60]},...
%                       'go',{'box','fancy'})
%
%    See also: PLOT0, PLOT1, PLOT2, RECORDSECTION, MAPEVENTGRID

%     Version History:
%        Dec.  2, 2009 - initial version
%        Dec.  8, 2009 - event grid plotting now MAPSTATIONS2
%        Mar.  1, 2010 - update for new checking state function names
%        May   7, 2010 - changed name to MAPSTATIONS
%        June 21, 2010 - use gshhs coastline & boarders
%        June 23, 2010 - proj & proj opt args
%        June 26, 2010 - opt/val pair inputs, several more options,
%                        plots all stations (undefined are set to NaN)
%        July 22, 2010 - uses MAPLOCATIONS
%        Aug. 21, 2010 - update undef usage
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 21, 2010 at 22:00 GMT

% todo:
% - msg on click features would be nice
%   - kname, lat, lon, depth, elevation
%   - a way to scatter dense on click (like googleearth)

% check nargin
error(nargchk(1,inf,nargin));
if(mod(nargin-1,2))
    error('seizmo:mapstations:badNumInputs',...
        'Unpaired Option/Value!');
end

% check data structure
error(seizmocheck(data));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt rest
try
    % get header info
    [stla,stlo,evla,evlo]=getheader(data,'stla','stlo','evla','evlo');
    
    % remove undefined
    badst=isnan(stla) | isnan(stlo);
    if(any(badst))
        warning('seizmo:mapstations:badLocation',...
            ['Station location not set for Records:\n' ...
            sprintf('%d ',find(badst))]);
    end

    % plot locations
    ax=maplocations('st',[stla stlo],'ev',[evla evlo],varargin{:});
    
    % output
    if(nargout); varargout{1}=ax; end

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
