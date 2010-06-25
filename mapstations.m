function [varargout]=mapstations(data,varargin)
%MAPSTATIONS    Map station/earthquake locations of SEIZMO records
%
%    Usage:    mapstations(data)
%              mapstations(data,gshhs)
%              mapstations(data,gshhs,proj)
%              mapstations(data,gshhs,proj,'opt1',val1,...,'optN',valN)
%              mapstations(...,ax)
%              ax=mapstations(...)
%
%    Description: MAPSTATIONS(DATA) creates a map showing the station and
%     earthquake locations stored in the headers of records of SEIZMO
%     struct DATA.  The map is a global map using the Hammer projection.
%     Stations are plotted as yellow circles and events are plotted as
%     5-pointed stars.
%
%     MAPSTATIONS(DATA,GSHHS) sets the GSHHS coastline and political
%     boundaries resolution.  The values can be 'c', 'l', 'i', 'h', 'f'.
%     The default GSHHS resolution is 'l' (low).
%
%     MAPSTATIONS(DATA,GSHHS,PROJ) defines the map projection.  See
%     M_PROJ('SET') for possible projections.  The default PROJ is
%     'Hammer-Aitoff'.
%
%     MAPSTATIONS(DATA,GSHHS,PROJ,'OPT1',VAL1,...,'OPTN',VALN) passes
%     additional options to M_PROJ.  See M_PROJ('get',PROJ) for a list of
%     possible options for this projection.  The default is no options,
%     which will create a global map.
%     
%     MAPSTATIONS(...,AX) uses the axes given by handle AX for the map.
%     This is the last input.  The default draws the map in a new figure.
%
%     AX=MAPSTATIONS(DATA) returns the axes handle for the map.
%
%    Notes:
%
%    Examples:
%     Show locations of stations in a dataset:
%      mapstations(data);
%
%    See also: PLOT0, PLOT1, PLOT2, RECORDSECTION, MAPSTATIONS2,
%              MAPCLUSTERS

%     Version History:
%        Dec.  2, 2009 - initial version
%        Dec.  8, 2009 - event grid plotting now MAPSTATIONS2
%        Mar.  1, 2010 - update for new checking state function names
%        May   7, 2010 - changed name to MAPSTATIONS
%        June 21, 2010 - use gshhs coastline & boarders
%        June 23, 2010 - proj & proj opt args
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 23, 2010 at 18:50 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% access to m_map globals for map boundaries
global MAP_VAR_LIST

% ocean/land colors
ocean=[0.3 0.6 1];
land=[0.4 0.6 0.2];
border=[0.5 0 0];

% check data structure
[h,idx]=versioninfo(data);

% get undefined values
undef=getsubfield(h,'undef','ntype').';
undef=undef(idx);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt rest
try
    % handle additional inputs
    skip=false;
    if(numel(varargin))
        % test if resolution or handle
        if(ischar(varargin{1}) || isempty(varargin{1}))
            % resolution
            if(isempty(varargin{1}))
                gshhs='l';
            else
                gshhs=varargin{1};
            end
            varargin(1)=[];
        else
            % default map
            gshhs='c';
            proj='hammer';
            varargin={'clon' 0};
            ax=varargin{1};
            varargin(1)=[];
            skip=true;
            
            % check if really a handle, otherwise ignore
            if(isempty(ax) || ~isscalar(ax) || ~isreal(ax) ...
                    || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
                % new figure
                figure;
                ax=gca;
            else
                axes(ax);
            end
            
            % if any more inputs then error
            if(numel(varargin))
                error('seizmo:mapstations:badInputs',...
                    'H must be the final input arg!');
            end
        end
    else
        % new figure and default map
        figure;
        ax=gca;
        gshhs='c';
        proj='hammer';
        varargin={'clon' 0};
        skip=true;
    end
    if(~skip && numel(varargin))
        % test if projection or handle
        if(ischar(varargin{1}) || isempty(varargin{1}))
            % projection
            if(isempty(varargin{1}))
                proj='hammer';
            else
                proj=varargin{1};
            end
            varargin(1)=[];
            
            % extract axes handle if it is there
            if(mod(numel(varargin),2))
                ax=varargin{end};
                varargin(end)=[];
                
                % check if really a handle, otherwise ignore
                if(isempty(ax) || ~isscalar(ax) || ~isreal(ax) ...
                        || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
                    % new figure
                    figure;
                    ax=gca;
                else
                    axes(ax);
                end
            else
                % new figure
                figure;
                ax=gca;
            end
        else % handle
            % default map
            proj='hammer';
            varargin={'clon' 0};
            ax=varargin{1};
            varargin(1)=[];
            
            % check if really a handle, otherwise ignore
            if(isempty(ax) || ~isscalar(ax) || ~isreal(ax) ...
                    || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
                % new figure
                figure;
                ax=gca;
            else
                axes(ax);
            end
            
            % if any more inputs then error
            if(numel(varargin))
                error('seizmo:mapstations:badInputs',...
                    'H must be the final input arg!');
            end
        end
    elseif(~skip)
        % new figure and default map
        figure;
        ax=gca;
        proj='hammer';
        varargin={'clon' 0};
    end
    
    % plot map
    m_proj(proj,varargin{:});
    set(ax,'color',ocean);
    m_gshhs([gshhs 'c'],'patch',land);
    m_gshhs([gshhs 'b'],'color',border);
    m_grid();
    
    % hackery to color oceans at large when the above fails
    set(findobj(ax,'tag','m_grid_color'),'facecolor',ocean);

    % get header info
    [stla,stlo,evla,evlo]=getheader(data,'stla','stlo','evla','evlo');
    
    % remove undefined
    badst=stla==undef | stlo==undef;
    stla(badst)=[]; stlo(badst)=[];
    if(any(badst))
        warning('seizmo:mapstations:badLocation',...
            ['Station location not set for Records:\n' ...
            sprintf('%d ',find(badst))]);
    end
    badev=evla==undef | evlo==undef;
    evla(badev)=[]; evlo(badev)=[];

    % get unique locations
    stlalo=unique([stla stlo],'rows');
    evlalo=unique([evla evlo],'rows');

    % wrap longitudes to plot
    while(any(abs(stlalo(:,2)-mean(MAP_VAR_LIST.longs))>180))
        stlalo(stlalo(:,2)<MAP_VAR_LIST.longs(1),2)=...
            stlalo(stlalo(:,2)<MAP_VAR_LIST.longs(1),2)+360;
        stlalo(stlalo(:,2)>MAP_VAR_LIST.longs(2),2)=...
            stlalo(stlalo(:,2)>MAP_VAR_LIST.longs(2),2)-360;
    end
    while(any(abs(evlalo(:,2)-mean(MAP_VAR_LIST.longs))>180))
        evlalo(evlalo(:,2)<MAP_VAR_LIST.longs(1),2)=...
            evlalo(evlalo(:,2)<MAP_VAR_LIST.longs(1),2)+360;
        evlalo(evlalo(:,2)>MAP_VAR_LIST.longs(2),2)=...
            evlalo(evlalo(:,2)>MAP_VAR_LIST.longs(2),2)-360;
    end

    % plot locations
    axes(ax);
    hold on
    h=m_scatter(evlalo(:,2),evlalo(:,1),200,'r','filled','p',...
        'markeredgecolor','k');
    set(h,'tag','events');
    h=m_scatter(stlalo(:,2),stlalo(:,1),[],'y','filled',...
        'markeredgecolor','k');
    set(h,'tag','stations');
    hold off
    
    % return figure handle
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
