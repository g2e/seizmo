function [varargout]=plotstations(data,fh)
%PLOTSTATIONS    Plots station and earthquake locations of SEIZMO records
%
%    Usage:    plotstations(data)
%              plotstations(data,h)
%              h=plotstations(...)
%
%    Description: PLOTSTATIONS(DATA) creates a map showing the station and
%     earthquake locations stored in the headers of records of SEIZMO
%     struct DATA.  The map is a global map using the Hammer projection.
%     Stations are plotted as yellow circles and events are plotted as
%     5-pointed stars.
%     
%     PLOTSTATIONS(DATA,H) uses the figure handle H for the map.
%
%     H=PLOTSTATIONS(DATA) returns the figure handle of the map figure.
%
%    Notes:
%
%    Examples:
%     Show locations of stations in a dataset:
%      h=plotstations(data);
%
%    See also: PLOT0, PLOT1, PLOT2, RECORDSECTION, PLOTSTATIONS2

%     Version History:
%        Dec.  2, 2009 - initial version
%        Dec.  8, 2009 - event grid plotting now PLOTSTATIONS2
%        Mar.  1, 2010 - update for new checking state function names
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2010 at 01:45 GMT

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
    % get figure
    if(nargin==2 && isscalar(fh) && fh==fix(fh))
        figure(fh);
    else
        fh=figure;
    end
    
    % plot map
    clon=0; % the default center longitude
    alon=[clon-180 clon+180]; % longitude boundaries
    m_proj('hammer','clon',clon);
    %m_proj('get')
    set(gca,'color',[0.6 0.96 1]);
    m_grid('xticklabels',[],'ytick',-90:15:90,'xtick',-180:15:180);
    m_coast('patch',[0.6 1 0.6]);

    % hackery to color oceans at large
    child=get(gca,'children');
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
    m_line(stlalo(:,2),stlalo(:,1),'marker','o',...
       'markerfacecolor','y','markersize',6,'markeredgecolor','k',...
       'linestyle','none');
    m_line(evlalo(:,2),evlalo(:,1),'marker','p',...
       'markerfacecolor','r','markersize',14,'markeredgecolor','k',...
       'linestyle','none');
    
    % return figure handle
    if(nargout); varargout{1}=fh; end

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

end
