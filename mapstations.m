function [varargout]=mapstations(data,han)
%MAPSTATIONS    Map station/earthquake locations of SEIZMO records
%
%    Usage:    mapstations(data)
%              mapstations(data,h)
%              h=mapstations(...)
%
%    Description: MAPSTATIONS(DATA) creates a map showing the station and
%     earthquake locations stored in the headers of records of SEIZMO
%     struct DATA.  The map is a global map using the Hammer projection.
%     Stations are plotted as yellow circles and events are plotted as
%     5-pointed stars.
%     
%     MAPSTATIONS(DATA,H) uses the axes given by handle H for the map.
%
%     H=MAPSTATIONS(DATA) returns the axes handle for the map.
%
%    Notes:
%
%    Examples:
%     Show locations of stations in a dataset:
%      h=mapstations(data);
%
%    See also: PLOT0, PLOT1, PLOT2, RECORDSECTION, MAPSTATIONS2,
%              MAPCLUSTERS

%     Version History:
%        Dec.  2, 2009 - initial version
%        Dec.  8, 2009 - event grid plotting now MAPSTATIONS2
%        Mar.  1, 2010 - update for new checking state function names
%        May   7, 2010 - changed name to MAPSTATIONS
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May   7, 2010 at 17:00 GMT

% todo:

% check nargin
msg=nargchk(1,2,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
[h,idx]=versioninfo(data);

% get undefined values
undef=getsubfield(h,'undef','ntype').';
undef=undef(idx);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt rest
try
    % check axes handle
    if(nargin==2 && ~isempty(han) && (~isreal(han) || ~isscalar(han)))
        error('seizmo:mapclusters:badInput',...
            'H must be an axis handle!');
    end
    
    % make new figure if no axes handle passed
    if(nargin<2 || isempty(han) || ~ishandle(han))
        figure();
        han=gca;
    else
        axes(han);
    end
    
    % plot map
    clon=0; % the default center longitude
    alon=[clon-180 clon+180]; % longitude boundaries
    m_proj('hammer','clon',clon);
    title('STATION MAP','fontweight','bold')
    %m_proj('get')
    set(han,'color',[0.6 0.96 1]);
    m_grid('xticklabels',[],'ytick',-90:15:90,'xtick',-180:15:180);
    m_coast('patch',[0.6 1 0.6]);

    % hackery to color oceans at large
    child=get(han,'children');
    try
        set(child(end),'facecolor',[0.6 0.96 1]);
    catch
        set(child(end-1),'facecolor',[0.6 0.96 1]);
    end

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
    while(any(abs(stlalo(:,2)-clon)>180))
        stlalo(stlalo(:,2)<alon(1),2)=stlalo(stlalo(:,2)<alon(1),2)+360;
        stlalo(stlalo(:,2)>alon(2),2)=stlalo(stlalo(:,2)>alon(2),2)-360;
    end
    while(any(abs(evlalo(:,2)-clon)>180))
        evlalo(evlalo(:,2)<alon(1),2)=evlalo(evlalo(:,2)<alon(1),2)+360;
        evlalo(evlalo(:,2)>alon(2),2)=evlalo(evlalo(:,2)>alon(2),2)-360;
    end

    % plot locations
    axes(han);
    hold on
    m_scatter(stlalo(:,2),stlalo(:,1),[],'y','filled',...
        'markeredgecolor','k');
    m_scatter(evlalo(:,2),evlalo(:,1),200,'r','filled','p',...
        'markeredgecolor','k');
    hold off
    
    % return figure handle
    if(nargout); varargout{1}=han; end

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

end
