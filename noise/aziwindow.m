function [varargout]=aziwindow(data,azwin,lagrng,degrng)
%AZIWINDOW    Azimuthal window record section of correlograms
%
%    Usage:    aziwindow(data,aziwin)
%              aziwindow(data,aziwin,lagrng,degrng)
%
%    Description:
%     
%     AZISWINDOW(DATA,AZIWIN) plots correlograms as a
%     function of
%     distance for a series of azimuthal
%     windows.  The default azimuthal window is 10 degrees wide and
%     increments at 1 degree.  The record section plot extends from -500 to
%     500 seconds in lag time and the distance range is -1deg to 11 deg.
%     To adjust these values see the remaining usage forms.  This also
%     plots a station map and a "look" window.  The map plotting can be
%     slow so expect to wait 10 or so minutes for this to finish.
%
%     MOV=AZIWINDOW(DATA,AZIWIN,LAGRNG,DEGRNG) adjusts the record section
%     plot limits.  Default LAGRNG is [-500 500] and DEGRNG is [-1 11].
%
%    Notes:
%
%    Examples:
%
%    See also: AZISWEEP, CORRELATE

%     Version History:
%        Sep.  7, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  7, 2010 at 14:15 GMT

% todo:

% check nargin
error(nargchk(1,4,nargin));

% check dataset
versioninfo(data,'dep');

% default/check optional inputs
if(nargin<2 || isempty(azwin)); azwin=[0 360]; end
if(nargin<3 || isempty(lagrng)); lagrng=[-500 500]; end
if(nargin<4 || isempty(degrng)); degrng=[-1 11]; end
if(numel(azwin)~=2 || ~isreal(azwin) || azwin(1)>azwin(2))
    error('seizmo:aziwindow:badInput',...
        'AZIWIN must be [AZMIN AZMAX] in degrees!');
end

% get necessary info
az=getheader(data,'az');

% find stations in window
in=(az>=azwin(1) & az<=azwin(2)) ...
    | (az>=(azwin(1)+360) & az<=(azwin(2)+360)) ...
    | (az>=(azwin(1)-360) & az<=(azwin(2)-360));

% skip if no records
if(~sum(in))
    error('seizmo:aziwindow:badInput',...
        'No pair in azimuthal window!');
end

% make plot
[ax1,ax2]=plot_noise_recsec(data(in),lagrng,degrng,azwin);

% output
if(nargout); varargout={ax1 ax2}; end

end

function [rs_ax,map_ax]=plot_noise_recsec(...
    data,lagrng,degrng,azrng)

% record section plot
rs_ax=recordsection(data,...
    'title',{['Azimuthal Window: ' num2str(azrng(1)) ' to ' ...
    num2str(azrng(2)) ' degrees'] ...
    [num2str(numel(data)) ' Station Pairs']},...
    'colormap','k',...
    'linewidth',1,...
    'xlim',lagrng,...
    'ylim',degrng,...
    'fgcolor','k','bgcolor','w',...
    'fontsize',10,...
    'ylabel','Distance (deg)',...
    'xlabel','Time (s)');
set(rs_ax,'linewidth',2);
grid(rs_ax,'off');

% stereographic station map
[st,ev]=getheader(data,'st','ev');
loc=unique([st; ev],'rows');
[clat,clon]=arraycenter(loc(:,1),loc(:,2));
%dist=max(sphericalinv(clat,clon,loc(:,1),loc(:,2)));
pos=get(rs_ax,'position');
map_ax=axes('position',[pos(1:2)+[0 .05] 0.3 0.3],...
    'parent',get(rs_ax,'parent'));
mapstations(data,'proj','stereo',...
    'po',{'lat',clat,'lon',clon,'rad',1.1*max(degrng)},...
    'go',{'xticklabels',[],'yticklabels',[]},...
    'gshhs','c',...
    'fgcolor','k',...
    'evms',75,...
    'axis',map_ax);

% add lines between station pairs
[lat,lon]=gcarc2latlon(ev(:,1),ev(:,2),st(:,1),st(:,2),20);
m_line(lon',lat','color','k')

% bold fonts
setfonts([rs_ax map_ax],'fontweight','bold');

end
