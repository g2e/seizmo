function [varargout]=azisweep(data,lagrng,degrng,azwin,azstep)
%AZISWEEP    Sliding azimuthal window record section of correlograms
%
%    Usage:    mov=azisweep(data)
%              mov=azisweep(data,lagrng,degrng)
%              mov=azisweep(data,lagrng,degrng,azwin,azstep)
%              mov=azisweep(data,lagrng,degrng,azwin,azstep,ax)
%
%    Description:
%     MOV=AZISWEEP(DATA) creates a Matlab movie MOV by drawing correlograms
%     as a function of distance for a series of azimuthal windows.  The
%     default azimuthal window is 10 degrees wide and steps at 1 degree.
%     The record section plot extends from -500 to 500 seconds in lag time
%     and the distance range is -1 to 11 degrees.  To adjust these values
%     see the remaining usage forms.  Also included is a small map inset on
%     the lower left showing the station pairs in the main window.  The
%     source stations are drawn as red stars while the receiver stations
%     are yellow circles.  Lines connecting the 2 stations designates a
%     pair.  The map plotting can be terribly slow so expect to wait 10 or
%     so minutes for this to finish.
%
%     MOV=AZISWEEP(DATA,LAGRNG,DEGRNG) adjusts the record section plot
%     limits.  Default LAGRNG is [-500 500] and DEGRNG is [-1 11].  Records
%     outside these limits but withing the azimuthal range are still drawn
%     in the map.
%
%     MOV=AZISWEEP(DATA,LAGRNG,DEGRNG,AZWIN,AZSTEP) adjusts the azimuth
%     window properties.  The default AZWIN is 10deg and AZSTEP is 1deg.
%
%     MOV=AZISWEEP(DATA,LAGRNG,DEGRNG,AZWIN,AZSTEP,AX) sets the axis to
%     draw the plot in to AX.  If AX contains 2 handles then the second
%     axis is used to draw the station map.  Note that both axes must share
%     the same parent (the second handle is ignored if they do not).
%
%    Notes:
%
%    Examples:
%     % Speed things up by increasing the azstep:
%     mov=azisweep(data,[],[],15,5);
%
%    See also: AZIWINDOW, CORRELATE, NOISE

%     Version History:
%        June 24, 2010 - initial version
%        Sep. 16, 2010 - doc update, fixed plotting (much faster now),
%                        fixed bug in azimuthal windowing and labeling
%        Feb. 11, 2011 - minor mlint fix
%        Apr.  3, 2012 - minor doc update, use seizmocheck
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 14:15 GMT

% todo:

% check nargin
error(nargchk(1,6,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% default/check optional inputs
if(nargin<2 || isempty(lagrng)); lagrng=[-500 500]; end
if(nargin<3 || isempty(degrng)); degrng=[-1 11]; end
if(nargin<4 || isempty(azwin)); azwin=10; end
if(nargin<5 || isempty(azstep)); azstep=1; end
if(nargin<6); ax=[]; end
if(numel(lagrng)~=2 || ~isreal(lagrng) || lagrng(1)>lagrng(2))
    error('seizmo:azisweep:badInput',...
        'LAGRNG must be [LAGMIN LAGMAX] in seconds!');
end
if(numel(degrng)~=2 || ~isreal(degrng) || degrng(1)>degrng(2))
    error('seizmo:azisweep:badInput',...
        'DEGRNG must be [DEGMIN DEGMAX] in degrees!');
end
if(~isscalar(azwin) || ~isreal(azwin) || azwin<=0)
    error('seizmo:azisweep:badInput',...
        'AZIWIN must be a positive real scalar (in degrees!)');
end
if(~isscalar(azstep) || ~isreal(azstep) || azstep<=0)
    error('seizmo:azisweep:badInput',...
        'AZISTEP must be a positive real scalar (in degrees!)');
end
if(numel(ax)>2 || ~isreal(ax))
    error('seizmo:aziwindow:badInput',...
        'AX must be vector of 1 or 2 handles!');
end

% do we need to save the movie frames?
makemovie=false;
if(nargout); makemovie=true; end

% get necessary info
[az,st,ev,dep]=getheader(data,'az','st','ev','dep');

% determine radius of the station map
loc=unique([st; ev],'rows');
[clat,clon]=arraycenter(loc(:,1),loc(:,2));
r=max(sphericalinv(clat,clon,loc(:,1),loc(:,2)));
r=1.25*r; % go a bit farther than the most distant station

% get maximum amplitude so we can scale azimuthal windows
% to preserve the relative amplitudes across frames
normmax=1/3; % this is what plot_noise_recsec uses
ampmax=max(abs(dep(:)));

% loop over azimuths
frame=1;
for k=0:azstep:360
    % find stations in window
    in=(az>=(k-azwin/2) & az<=(k+azwin/2)) ...
        | (az>=(k-azwin/2+360) & az<=(k+azwin/2+360)) ...
        | (az>=(k-azwin/2-360) & az<=(k+azwin/2-360));
    
    % skip if no records
    if(~sum(in)); continue; end
    
    % get normalizer for this window
    tmpdep=dep(in,:);
    newnormmax=normmax*max(abs(tmpdep(:)))/ampmax;
    
    % make plot
    if(frame==1)
        ax=plot_noise_recsec(data(in),...
            [k-azwin/2 k+azwin/2],lagrng,degrng,newnormmax,...
            ev(in,:),st(in,:),clat,clon,r,ax);
    else
        % redraw the records, replace points/lines on the map
        update_noise_recsec(data(in),...
            [k-azwin/2 k+azwin/2],lagrng,degrng,newnormmax,...
            ev(in,:),st(in,:),ax);
    end
    
    % save frame
    if(makemovie); varargout{1}(frame)=getframe(get(ax(1),'parent')); end
    frame=frame+1;
end

end

function [ax]=plot_noise_recsec(data,...
    azrng,lagrng,degrng,normmax,ev,st,clat,clon,r,ax)

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
    'normmax',normmax,...
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
        || ~strcmp('axes',get(ax(2),'type')) ...
        || get(ax(1),'parent')~=get(ax(2),'parent'))
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
mapstations(data,...
    'proj','stereo',...
    'po',{'lat',clat,'lon',clon,'rad',r},...
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


function []=update_noise_recsec(data,azrng,lagrng,degrng,normmax,ev,st,ax)

% crash if axes no longer exist
if(any(~ishandle(ax)))
    error('seizmo:azisweep:deletedHandle',...
        'Axis was deleted!  Can not continue making movie.');
end

% redraw record section plot
ax(1)=recordsection(data,...
    'title',{['Azimuthal Window: ' num2str(azrng(1)) '^o to ' ...
    num2str(azrng(2)) '^o'] ...
    [num2str(numel(data)) ' Station Pairs']},...
    'colormap','k',...
    'linewidth',0.1,...
    'xlim',lagrng,...
    'ylim',degrng,...
    'normmax',normmax,...
    'fgcolor','k',...
    'bgcolor','w',...
    'fontsize',10,...
    'ylabel','Station Pair Distance (^o)',...
    'xlabel','Lag Time (s)',...
    'ax',ax(1));
set(ax(1),'linewidth',2);
grid(ax(1),'off');

% station and line handles
axes(ax(2));
he=findobj(ax(2),'tag','events');
hs=findobj(ax(2),'tag','stations');
hl=findobj(ax(2),'tag','stpairlines');

% delete lines
delete(hl);

% access to m_map globals for map boundaries
global MAP_VAR_LIST

% wrap longitudes to plot
while(any(abs(st(:,2)-mean(MAP_VAR_LIST.longs))>180))
    st(st(:,2)<MAP_VAR_LIST.longs(1),2)=...
        st(st<MAP_VAR_LIST.longs(1),2)+360;
    st(st(:,2)>MAP_VAR_LIST.longs(2),2)=...
        st(st(:,2)>MAP_VAR_LIST.longs(2),2)-360;
end
while(any(abs(ev(:,2)-mean(MAP_VAR_LIST.longs))>180))
    ev(ev(:,2)<MAP_VAR_LIST.longs(1),2)=...
        ev(ev(:,2)<MAP_VAR_LIST.longs(1),2)+360;
    ev(ev(:,2)>MAP_VAR_LIST.longs(2),2)=...
        ev(ev(:,2)>MAP_VAR_LIST.longs(2),2)-360;
end

% replace station info
[x,y]=m_ll2xy(ev(:,2),ev(:,1));
set(he,'xdata',x,'ydata',y);
[x,y]=m_ll2xy(st(:,2),st(:,1));
set(hs,'xdata',x,'ydata',y);

% draw lines
[lat,lon]=gcarc2latlon(ev(:,1),ev(:,2),st(:,1),st(:,2),20);
hl=m_line(lon',lat','color','k','linewidth',0.1);
set(hl,'tag','stpairlines');

% bold fonts
setfonts(ax,'fontweight','bold');

end


