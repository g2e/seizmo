function [varargout]=mapstations(data,varargin)
%MAPSTATIONS    Map station/earthquake locations of SEIZMO records
%
%    Usage:    mapstations(data)
%              mapstations(...,'mmap_opt1',mmap_val1,...)
%              ax=mapstations(...)
%
%    Description:
%     MAPSTATIONS(DATA) creates a map showing the station and earthquake
%     locations stored in the headers of records of SEIZMO struct DATA
%     using MMAP defaults.  So stations are plotted as yellow circles and
%     events are plotted as 5-pointed red stars.
%
%     MAPSTATIONS(...,'MMAP_OPT1',MMAP_VAL1,...) passes additional options
%     on to MMAP to alter the map/symbols.
%
%     AX=MAPSTATIONS(...) returns the axes drawn in.
%
%    Notes:
%     - Station symbols are tagged as 'stations'.
%     - Event symbols are tagged as 'events'.
%
%    Examples:
%     % Show locations of stations in a dataset:
%     mapstations(data);
%
%     % This replaces the old MAPSTATIONS2:
%     ax=mapstations(data);
%     ev=getheader(data(1),'ev');
%     mapeventgrid(ax,ev(1),ev(2));
%
%     % This replaces the old MAPCLUSTERS:
%     ax=mapstations(data);
%     h=findobj(ax,'tag','stations');
%     set(h,'cdata',grp.color(grp.T,:));
%
%     % Show Africa stations in a map with fancy border:
%     mapstations(data,'po',{'lat',[-40 40],'lon',[-30 60]},...
%                      'go',{'box','fancy'})
%
%    See also: PLOT0, PLOT1, PLOT2, RECORDSECTION, MAPEVENTGRID, MMAP

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
%        Feb. 10, 2011 - update for maplocations=>mmap
%        Apr.  3, 2012 - minor doc update
%        Aug. 28, 2013 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 28, 2013 at 22:00 GMT

% todo:

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
    ax=mmap('st',[stla stlo],'ev',[evla evlo],varargin{:});
    
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
