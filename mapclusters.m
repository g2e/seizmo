function [varargout]=mapclusters(data,grp,han)
%MAPCLUSTERS    Map station locations colored by cluster
%
%    Usage:    mapclusters(data,grp)
%              mapclusters(data,grp,h)
%              h=mapclusters(...)
%
%    Description: MAPCLUSTERS(DATA,GRP) plots the station locations for
%     records in SEIZMO struct DATA with the group coloring stored in the
%     struct GRP output from USERCLUSTER.  Useful for visualizing if
%     clustering of waveforms is regionally constrained.  The plot is drawn
%     in a new figure.
%
%     MAPCLUSTERS(DATA,GRP,H) makes the cluster map in the axes given by
%     handle H.  If H is empty, a new figure is created for the plot.
%
%     H=MAPCLUSTERS(...) returns the axes handle for the cluster map.
%
%    Notes:
%
%    Examples:
%     The typical ordering (note USERCLUSTER can call MAPCLUSTERS too):
%      grp=usercluster(data);
%      mapclusters(data,grp);
%
%    See also: MAPSTATIONS, USERCLUSTER

%     Version History:
%        May   7, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May   7, 2010 at 16:45 GMT

% todo:

% check nargin
msg=nargchk(2,3,nargin);
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
    % check inputs
    if(~isstruct(grp) || ~all(ismember({'T' 'color'},fieldnames(grp))))
        error('seizmo:mapclusters:badInput',...
            'GRP is not formatted correctly!');
    elseif(nargin==3 && ~isempty(han) && (~isreal(han) || ~isscalar(han)))
        error('seizmo:mapclusters:badInput',...
            'H must be an axis handle!');
    end
    
    % make new figure if no axes handle passed
    if(nargin<3 || isempty(han) || ~ishandle(han))
        figure();
        han=gca;
    else
        axes(han);
    end
    
    % plot map
    clon=0; % the default center longitude
    alon=[clon-180 clon+180]; % longitude boundaries
    m_proj('hammer','clon',clon);
    title('CLUSTER MAP','fontweight','bold')
    %m_proj('get')
    set(han,'color',[0.6 0.96 1]); % oceans
    m_grid('xticklabels',[],'ytick',-90:15:90,'xtick',-180:15:180);
    m_coast('patch',[0.6 1 0.6]); % continents

    % hackery to completely color oceans
    child=get(han,'children');
    try
        set(child(end),'facecolor',[0.6 0.96 1]);
    catch
        set(child(end-1),'facecolor',[0.6 0.96 1]);
    end

    % get header info
    [stla,stlo,evla,evlo]=getheader(data,'stla','stlo','evla','evlo');
    
    % get colors
    stcolor=grp.color(grp.T,:);
    
    % remove undefined
    badst=stla==undef | stlo==undef;
    stla(badst)=[]; stlo(badst)=[]; stcolor(badst,:)=[];
    if(any(badst))
        warning('seizmo:mapclusters:badLocation',...
            ['Station location not set for Records:\n' ...
            sprintf('%d ',find(badst))]);
    end
    badev=evla==undef | evlo==undef;
    evla(badev)=[]; evlo(badev)=[];

    % combine lat/lon
    stlalo=[stla stlo];
    evlalo=[evla evlo];

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
    m_scatter(stlalo(:,2),stlalo(:,1),[],stcolor,'filled',...
        'markeredgecolor','k');
    m_scatter(evlalo(:,2),evlalo(:,1),200,'r','filled','p',...
        'markeredgecolor','k');
    hold off
    
    % return axis handle
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
