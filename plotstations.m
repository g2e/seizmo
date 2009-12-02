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
%    See also: PLOT0, PLOT1, PLOT2, RECORDSECTION

%     Version History:
%        Dec.  2, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec.  2, 2009 at 22:45 GMT

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
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

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
    m_grid('xticklabels',[]);
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

    % event-based gridding parameters
    %ranges=90:10:170;
    %azims=0:30:360;

    % plot azimuthal lines
    %[azisla,azislo]=sphericalfwd(evlalo(1,1),evlalo(1,2),ranges(1),azims);
    %[aziela,azielo]=sphericalfwd(evlalo(1,1),evlalo(1,2),ranges(end),azims);
    %[azimla,azimlo]=...
    %    sphericalfwd(evlalo(1,1),evlalo(1,2),sum(ranges([1 end]))/2,azims);
    %[lat,lon]=gcarc2latlon(azisla,azislo,azimla,azimlo,200);
    %m_line(lon',lat','color','b','linewi',2)
    %[lat,lon]=gcarc2latlon(aziela,azielo,azimla,azimlo,200);
    %m_line(lon',lat','color','b','linewi',2)

    % plot range rings
    %m_range_ring(evlalo(1,2),evlalo(1,1),6371*pi/180*ranges,200,...
    %    'color','b','linewi',2);
    
    % return figure handle
    if(nargout); varargout{1}=fh; end

    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

end
