function [varargout]=mapstations(data,varargin)
%MAPSTATIONS    Map station/earthquake locations of SEIZMO records
%
%    Usage:    mapstations(data)
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
%     struct DATA.  The map is a global map using the Hammer projection.
%     Stations are plotted as yellow circles and events are plotted as
%     5-pointed stars.
%
%     MAPSTATIONS(...,'GSHHS',RES,...) sets the GSHHS coastline and
%     political boundaries resolution.  The values can be 'c', 'l', 'i',
%     'h', 'f', or 'o' (for 'off').  The default GSHHS resolution is 'o'
%     (off) which calls M_COAST and does not draw political borders.
%
%     MAPSTATIONS(...,'PROJ',PROJ,...) defines the map projection.  See
%     M_PROJ('SET') for possible projections.  The default PROJ is
%     'Hammer-Aitoff'.
%
%     MAPSTATIONS(...,'PROJOPT',{'opt',val,...},...) passes additional
%     options to M_PROJ (like the lat/lon boundaries of the map).  The
%     options must be wrapped in a cell array!  See M_PROJ('get',PROJ) for
%     a list of possible options for the set projection (see 'PROJ' option
%     for the default projection and altering it).  The default is no
%     additional options, which will create a global map.
%
%     MAPSTATIONS(...,'GRIDOPT',{'opt',val,...},...) passes options to
%     M_GRID (like the lat/lon ticks of the map, etc).  The options must be
%     wrapped in a cell array!  See M_GRID('get') for a list of possible
%     options and M_GRID('set') for their defaults.  The default is no
%     options.
%
%     MAPSTATIONS(...,'FGCOLOR',COLOR,...) specifies the foreground color
%     of the map.  The default is 'w'.  If BGCOLOR is specified and FGCOLOR
%     is not, then FGCOLOR will be set using INVERTCOLOR.
%
%     MAPSTATIONS(...,'BGCOLOR',COLOR,...) specifies the background color
%     of the map.  The default is 'k'.  If FGCOLOR is specified and BGCOLOR
%     is not, then BGCOLOR will be set using INVERTCOLOR.
%
%     MAPSTATIONS(...,'SEA',COLOR,...) specifies the color of the sea in
%     the map.  The default is [.3 .6 1].
%
%     MAPSTATIONS(...,'LAND',COLOR,...) specifies the color of the land in
%     the map.  The default is [.4 .6 .2].
%
%     MAPSTATIONS(...,'BORDER',COLOR,...) specifies the color of the
%     political borders in the map.  The default is [.5 0 0].
%
%     MAPSTATIONS(...,'AXIS',AX,...)  sets the axes to draw in.  This is
%     useful for subplots, guis, etc.  The default draws the map in a new
%     figure.
%
%     AX=MAPSTATIONS(DATA) returns the axes handle for the map.
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 26, 2010 at 19:30 GMT

% todo:
% - msg on click features would be nice
%   - kname, lat, lon, depth, elevation
% - need a way to scatter dense on click (like googleearth)

% check nargin
error(nargchk(1,inf,nargin));
if(mod(nargin-1,2))
    error('seizmo:mapstations:badNumInputs',...
        'Unpaired Option/Value!');
end

% access to m_map globals for map boundaries
global MAP_VAR_LIST

% check data structure
[h,idx]=versioninfo(data);

% get undefined values
undef=getsubfield(h,'undef','ntype').';
undef=undef(idx);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt rest
try
    % option defaults
    varargin=[{'g' 'o' 'p' 'hammer' 'po' [] 'go' [] 'fg' [] 'bg' [] ...
        's' [.3 .6 1] 'l' [.4 .6 .2] 'b' [.5 0 0] 'a' []} varargin];
    
    % check options
    if(~iscellstr(varargin(1:2:end)))
        error('seizmo:mapstations:badOption',...
            'All Options must be specified with a string!');
    end
    for i=1:2:numel(varargin)
        % skip empty by default (but still checking option exists)
        skip=false;
        if(isempty(varargin{i+1})); skip=true; end
        val=varargin{i+1};
        
        % check option is available
        switch lower(varargin{i})
            case {'gshhs' 'res' 'g'}
                if(skip); continue; end
                if(ischar(val) && numel(val)==1 ...
                        && any(strcmpi(val,{'o' 'c' 'l' 'i' 'h' 'f'})))
                    gshhs=lower(val);
                else
                    error('seizmo:mapstations:badInput',...
                        'GSHHS option must be c, l, i, h, or f !');
                end
            case {'projection' 'proj' 'p'}
                if(skip); continue; end
                if(ischar(val) && ndims(val)==2 && size(val,1)==1)
                    proj=lower(val);
                else
                    error('seizmo:mapstations:badInput',...
                        'PROJECTION option must be a string!');
                end
            case {'projopt' 'popt' 'po'}
                if(skip)
                    popt={};
                elseif(iscell(val) && iscellstr(val(1:2:end)))
                    popt=val;
                else
                    error('seizmo:mapstations:badInput',...
                        ['PROJOPT option must be a cell array of ' ...
                        '''option''/value pairs!']);
                end
            case {'gridopt' 'gopt' 'go'}
                if(skip)
                    gopt={};
                elseif(iscell(val) && iscellstr(val(1:2:end)))
                    gopt=val;
                else
                    error('seizmo:mapstations:badInput',...
                        ['GRIDOPT option must be a cell array of ' ...
                        '''option''/value pairs!']);
                end
            case {'fgcolor' 'fg'}
                if(skip)
                    fg=[];
                elseif(ischar(val) ...
                        || (isreal(val) && isequal(size(val),[1 3])))
                    fg=val;
                else
                    error('seizmo:mapstations:badInput',...
                        'FGCOLOR must be a colorname or RGB triplet!');
                end
            case {'bgcolor' 'bg'}
                if(skip)
                    bg=[];
                elseif(ischar(val) ...
                        || (isreal(val) && isequal(size(val),[1 3])))
                    bg=val;
                else
                    error('seizmo:mapstations:badInput',...
                        'BGCOLOR must be a colorname or RGB triplet!');
                end
            case {'seacolor' 'sea' 's'}
                if(skip); continue; end
                if(ischar(val) ...
                        || (isreal(val) && isequal(size(val),[1 3])))
                    sea=val;
                else
                    error('seizmo:mapstations:badInput',...
                        'SEACOLOR must be a colorname or RGB triplet!');
                end
            case {'landcolor' 'land' 'l'}
                if(skip); continue; end
                if(ischar(val) ...
                        || (isreal(val) && isequal(size(val),[1 3])))
                    land=val;
                else
                    error('seizmo:mapstations:badInput',...
                        'LANDCOLOR must be a colorname or RGB triplet!');
                end
            case {'bordercolor' 'border' 'b'}
                if(skip); continue; end
                if(ischar(val) ...
                        || (isreal(val) && isequal(size(val),[1 3])))
                    border=val;
                else
                    error('seizmo:mapstations:badInput',...
                        'BORDERCOLOR must be a colorname or RGB triplet!');
                end
            case {'axis' 'ax' 'a'}
                if(skip)
                    ax=[];
                else
                    ax=val;
                end
            otherwise
                error('seizmo:mapstations:badOption',...
                    'Unknown Option: %s',varargin{i});
        end
    end
    
    % fix fg/bg colors
    if(isempty(fg))
        if(isempty(bg))
            fg='w'; bg='k';
        else
            fg=invertcolor(bg,true);
        end
    elseif(isempty(bg))
        bg=invertcolor(fg,true);
    end
    
    % convert colornames
    if(ischar(fg)); fg=name2rgb(fg); end
    if(ischar(bg)); bg=name2rgb(bg); end
    if(ischar(sea)); sea=name2rgb(sea); end
    if(ischar(land)); land=name2rgb(land); end
    if(ischar(border)); border=name2rgb(border); end
    
    % setup axis
    if(isempty(ax) || ~isscalar(ax) || ~isreal(ax) ...
            || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
        % new figure
        figure('color',bg);
        ax=gca;
    else
        axes(ax);
        h=get(ax,'children'); delete(h);
        h=findobj(get(get(ax,'parent'),'children'),'peer',ax); delete(h);
    end
    
    % plot map
    m_proj(proj,popt{:});
    set(ax,'color',sea);
    if(strcmpi(gshhs,'o'))
        m_coast('patch',land);
    else
        m_gshhs([gshhs 'c'],'patch',land);
        m_gshhs([gshhs 'b'],'color',border);
    end
    m_grid('color',fg,gopt{:});
    
    % hackery to color oceans at large when the above fails
    set(findobj(ax,'tag','m_grid_color'),'facecolor',sea);

    % get header info
    [stla,stlo,evla,evlo]=getheader(data,'stla','stlo','evla','evlo');
    
    % remove undefined
    badst=stla==undef | stlo==undef;
    stla(badst)=nan; stlo(badst)=nan;
    if(any(badst))
        warning('seizmo:mapstations:badLocation',...
            ['Station location not set for Records:\n' ...
            sprintf('%d ',find(badst))]);
    end
    badev=evla==undef | evlo==undef;
    evla(badev)=nan; evlo(badev)=nan;

    % wrap longitudes to plot
    while(any(abs(stlo-mean(MAP_VAR_LIST.longs))>180))
        stlo(stlo<MAP_VAR_LIST.longs(1))=...
            stlo(stlo<MAP_VAR_LIST.longs(1))+360;
        stlo(stlo>MAP_VAR_LIST.longs(2))=...
            stlo(stlo>MAP_VAR_LIST.longs(2))-360;
    end
    while(any(abs(evlo-mean(MAP_VAR_LIST.longs))>180))
        evlo(evlo<MAP_VAR_LIST.longs(1))=...
            evlo(evlo<MAP_VAR_LIST.longs(1))+360;
        evlo(evlo>MAP_VAR_LIST.longs(2))=...
            evlo(evlo>MAP_VAR_LIST.longs(2))-360;
    end

    % plot locations
    axes(ax);
    hold on
    h=m_scatter(evlo,evla,200,'r','filled','p','markeredgecolor','k');
    set(h,'tag','events');
    h=m_scatter(stlo,stla,[],'y','filled','markeredgecolor','k');
    set(h,'tag','stations');
    hold off
    
    % return figure handle
    set(ax,'tag','stationmap');
    if(nargout); varargout{1}=ax; end

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

end
