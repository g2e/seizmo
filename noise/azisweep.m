function [varargout]=azisweep(data,lagrng,degrng,azwin,azstep)
%AZISWEEP    Sliding azimuthal window record section of correlograms
%
%    Usage:    mov=azisweep(data)
%              mov=azisweep(data,lagrng,degrng,)
%              mov=azisweep(data,lagrng,degrng,azwin,azstep)
%
%    Description: MOV=AZISWEEP(DATA) creates a movie by plotting
%     correlograms as a function of distance for a series of azimuthal
%     windows.  The default azimuthal window is 10 degrees wide and
%     increments at 1 degree.  The record section plot extends from -500 to
%     500 seconds in lag time and the distance range is -1deg to 11 deg.
%     To adjust these values see the remaining usage forms.  This also
%     plots a station map and a "look" window.  The map plotting can be
%     slow so expect to wait 10 or so minutes for this to finish.
%
%     MOV=AZISWEEP(DATA,LAGRNG,DEGRNG) adjusts the record section plot
%     limits.  Default LAGRNG is [-500 500] and DEGRNG is [-1 11].
%
%     MOV=AZISWEEP(DATA,LAGRNG,DEGRNG,AZWIN,AZSTEP) adjusts the azimuth
%     window properties.  The default AZWIN is 10deg and AZSTEP is 1deg.
%
%    Notes:
%
%    Examples:
%     Speed things up by increasing the azstep:
%      mov=azisweep(data,[],[],15,5);
%
%    See also: CORRELATE, ROTATE_CORRELATIONS, REVERSE_CORRELATIONS

%     Version History:
%        June 24, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 24, 2010 at 14:15 GMT

% todo:

% check nargin
error(nargchk(1,5,nargin));

% check dataset
versioninfo(data,'dep');

% default/check optional inputs
if(nargin<2 || isempty(lagrng)); lagrng=[-500 500]; end
if(nargin<3 || isempty(degrng)); degrng=[-1 11]; end
if(nargin<4 || isempty(azwin)); azwin=10; end
if(nargin<5 || isempty(azstep)); azstep=1; end
if(~isscalar(azwin) || ~isreal(azwin) || azwin<=0)
    error('seizmo:azisweep:badInput',...
        'AZIWIN must be a positive real scalar (in degrees!)');
end
if(~isscalar(azstep) || ~isreal(azstep) || azstep<=0)
    error('seizmo:azisweep:badInput',...
        'AZISTEP must be a positive real scalar (in degrees!)');
end

% do we make the movie
makemovie=false;
if(nargout); makemovie=true; end

% get necessary info
az=getheader(data,'az');

% loop over azimuths
frame=1; fh=figure;
for k=0:azstep:360
    % find stations in window
    in=(az>(k-azwin/2) & az<(k+azwin/2)) ...
        | (az>(k-azwin/2+360) & az<(k-azwin/2+360)) ...
        | (az>(k-azwin/2-360) & az<(k-azwin/2-360));
    
    % skip if no records
    if(~sum(in)); continue; end
    
    % make plot
    if(frame==1)
        [ax1,ax2,ax3]=plot_noise_recsec(data(in),lagrng,degrng,...
            [k-azwin k+azwin],fh);
    else
        h=get(ax1,'children'); delete(h);
        h=get(ax2,'children'); delete(h);
        h=get(ax3,'children'); delete(h);
        delete(ax1); delete(ax2); delete(ax3);
        [ax1,ax2,ax3]=plot_noise_recsec(data(in),lagrng,degrng,...
            [k-azwin k+azwin],fh);
    end
    
    % save frame
    if(makemovie); varargout{1}(frame)=getframe(fh); frame=frame+1; end
end

end

function [rs_ax,az_ax,map_ax]=plot_noise_recsec(...
    data,lagrng,degrng,azrng,fh)

% record section plot
recordsection(data,...
    'fighandle',fh,...
    'title',{['Azimuthal Window: ' num2str(azrng(1)) '\circ to ' ...
    num2str(azrng(2)) '\circ'] [num2str(numel(data)) ' Station Pairs']},...
    'colormap','k',...
    'axislinewidth',2,...
    'recwidth',1,...
    'xlimits',lagrng,...
    'ylimits',degrng,...
    'fgcolor','k','bgcolor','w',...
    'fontweight','bold','fontsize',10,...
    'ylabel','Distance (\circ)',...
    'xlabel','Time (s)');
rs_ax=gca;

% azimuthal window plot
figure(fh);
theta=[azrng(1) azrng(1):diff(azrng)/50:azrng(2) azrng(2)]*pi/180;
rho=[0 ones(1,51) 0];
az_ax=axes('position',[0.7 0.03 0.3 0.3]);
h=mmpolar(theta,rho,...
    'style','compass',...
    'ttickvalue',0:30:330,...
    'rlimit',[0 1.05],...
    'tticklabelvisible','off',...
    'rticklabelvisible','off',...
    'grid','off');
hold on
patch(get(h,'xdata'), get(h,'ydata'),'b');
hold off

% stereographic station map
[st,ev]=getheader(data,'st','ev');
loc=unique([st; ev],'rows');
[clat,clon]=arraycenter(loc(:,1),loc(:,2));
dist=max(sphericalinv(clat,clon,loc(:,1),loc(:,2)));
figure(fh);
map_ax=axes('position',[0.03 0.03 0.3 0.3]);
mapstations(data,'proj','stereo',...
    'po',{'lat',clat,'lon',clon,'rad',dist*1.25},...
    'gshhs','c','axis',map_ax);

% add lines between station pairs
[lat,lon]=gcarc2latlon(ev(:,1),ev(:,2),st(:,1),st(:,2),20);
m_line(lon',lat','color','k')

end
