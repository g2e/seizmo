function [varargout]=aziwindow(data,azrng,lagrng,degrng,ax)
%AZIWINDOW    Azimuthal window record section of correlograms
%
%    Usage:    aziwindow(data,azirng)
%              aziwindow(data,azirng,lagrng,degrng)
%              aziwindow(data,azirng,lagrng,degrng,ax)
%              ax=aziwindow(...)
%
%    Description:
%     AZIWINDOW(DATA,AZIRNG) creates a plot showing the correlograms in
%     SEIZMO struct DATA ordered by station-pair distance.  The plot has
%     default limits of -1 to 11 in degree distance and +/-500s in time.
%     To adjust these values see the remaining usage forms.  The input
%     parameter AZIRNG is a 1x2 vector of [MINAZ MAXAZ] that limits the
%     data plotted to those station pairs with a source station to receiver
%     station azimuth within the range.  The default value is [0 360].
%     Also included is a small map inset on the lower left showing the
%     station pairs in the main window.  The source stations are drawn as
%     red stars while the receiver stations are yellow circles.  Lines
%     connecting the 2 stations designates a pair.
%
%     AZIWINDOW(DATA,AZIRNG,LAGRNG,DEGRNG) adjusts the record section
%     plot limits.  Default LAGRNG is [-500 500] and DEGRNG is [-1 11].
%
%     AZIWINDOW(DATA,AZIRNG, LAGRNG, DEGRNG,AX) sets the axis to draw the
%     plot in to AX.  This is useful for subplots, guis, etc.  If AX
%     contains 2 handles then the second axis is used to draw the station
%     map.
%
%     AX=AZIWINDOW(...) returns the axes of the record section and the map.
%
%    Notes:
%
%    Examples:
%     % Make two plots, one showing station North/South
%     % pairs while the other shows East/West pairs:
%     aziwindow(data,[-20 20])
%     aziwindow(data,[70 110])
%
%    See also: AZISWEEP, CORRELATE, NOISE

%     Version History:
%        Sep.  7, 2010 - initial version
%        Sep. 16, 2010 - major doc update, better checks and plotting
%        Apr.  3, 2012 - use seizmocheck
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 14:15 GMT

% todo:

% check nargin
error(nargchk(1,5,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% check header of dataset
data=checkheader(data,...
    'NONTIME_IFTYPE','ERROR',...
    'FALSE_LEVEN','ERROR',...
    'UNSET_ST_LATLON','ERROR',...
    'UNSET_EV_LATLON','ERROR',...
    'UNSET_DELAZ','ERROR');

% default optional inputs
if(nargin<2 || isempty(azrng)); azrng=[0 360]; end
if(nargin<3 || isempty(lagrng)); lagrng=[-500 500]; end
if(nargin<4 || isempty(degrng)); degrng=[-1 11]; end
if(nargin<5); ax=[]; end

% check inputs
if(numel(azrng)~=2 || ~isreal(azrng) || azrng(1)>azrng(2))
    error('seizmo:aziwindow:badInput',...
        'AZIRNG must be [AZMIN AZMAX] in degrees!');
end
if(numel(lagrng)~=2 || ~isreal(lagrng) || lagrng(1)>lagrng(2))
    error('seizmo:aziwindow:badInput',...
        'LAGRNG must be [LAGMIN LAGMAX] in seconds!');
end
if(numel(degrng)~=2 || ~isreal(degrng) || degrng(1)>degrng(2))
    error('seizmo:aziwindow:badInput',...
        'DEGRNG must be [DEGMIN DEGMAX] in degrees!');
end
if(numel(ax)>2 || ~isreal(ax))
    error('seizmo:aziwindow:badInput',...
        'AX must be vector of 1 or 2 handles!');
end

% get necessary info
az=getheader(data,'az');

% find stations in window
in=(az>=azrng(1) & az<=azrng(2)) ...
    | (az>=(azrng(1)+360) & az<=(azrng(2)+360)) ...
    | (az>=(azrng(1)-360) & az<=(azrng(2)-360));

% skip if no records
if(~sum(in))
    error('seizmo:aziwindow:badInput',...
        'No pair in azimuthal window!');
end

% make plot
ax=plot_noise_recsec(data(in),azrng,lagrng,degrng,ax);

% output
if(nargout); varargout={ax}; end

end

function [ax]=plot_noise_recsec(data,azrng,lagrng,degrng,ax)

% check first axis
if(isempty(ax) || ~isreal(ax) || ~ishandle(ax(1)) ...
        || ~strcmp('axes',get(ax(1),'type')))
    % new figure
    fh=figure('color','w');
    ax(1)=axes('parent',fh);
    nfflag=true;
else
    cla(ax(1),'reset');
    axes(ax(1));
    nfflag=false;
end

% record section plot
ax(1)=recordsection(data,...
    'title',{['Azimuthal Window: ' num2str(azrng(1)) '^o to ' ...
    num2str(azrng(2)) '^o'] ...
    [num2str(numel(data)) ' Station Pairs']},...
    'colormap','k',...
    'linewidth',0.1,...
    'xlim',lagrng,...
    'ylim',degrng,...
    'fgcolor','k',...
    'bgcolor','w',...
    'fontsize',10,...
    'ylabel','Station Pair Distance (^o)',...
    'xlabel','Lag Time (s)',...
    'ax',ax(1));
set(ax(1),'linewidth',2);
grid(ax(1),'off');
pos=get(ax(1),'position');

% check second axis
if(numel(ax)<2 || ~isreal(ax) || ~ishandle(ax(2)) ...
        || ~strcmp('axes',get(ax(2),'type')))
    % new axis
    if(nfflag)
        ax(2)=axes('parent',get(ax(1),'parent'),...
            'position',[pos(1:2)+[0 .05] 0.3 0.3]);
    else
        ax(2)=axes('parent',get(ax(1),'parent'),...
            'position',pos.*[1 1 .3 .3]+[0 1/15*pos(4) 0 0]);
    end
else
    % clear and bring to focus
    cla(ax(2),'reset');
    axes(ax(2));
end

% stereographic station map
[st,ev]=getheader(data,'st','ev');
loc=unique([st; ev],'rows');
[clat,clon]=arraycenter(loc(:,1),loc(:,2));
dist=max(sphericalinv(clat,clon,loc(:,1),loc(:,2)));
mapstations(data,...
    'proj','stereo',...
    'po',{'lat',clat,'lon',clon,'rad',1.25*dist},...
    'go',{'xticklabels',[],'yticklabels',[]},...
    'gshhs','i',...
    'fgcolor','k',...
    'evms',75,...
    'axis',ax(2));

% add lines between station pairs
[lat,lon]=gcarc2latlon(ev(:,1),ev(:,2),st(:,1),st(:,2),20);
hl=m_line(lon',lat','color','k','linewidth',0.1);
set(hl,'tag','stpairlines');

% bold fonts
setfonts(ax,'fontweight','bold');

end
