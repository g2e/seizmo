function [varargout]=plotgeofss(s,dblim,zerodb,varargin)
%PLOTGEOFSS    Plots frequency-slowness-position spectra
%
%    Usage:    plotgeofss(s)
%              plotgeofss(s,dblim)
%              plotgeofss(s,dblim,zerodb)
%              plotgeofss(s,dblim,zerodb,'param1',val1,...)
%
%    Description:
%     PLOTGEOFSS(S) plots the frequency-slowness-position power spectra
%     in struct S.  See a geofss function like GEOFSSXC for details on the
%     struct.  Note that the GEOFSS power spectra must be scalar in the
%     frequency and slowness dimensions (see GEOFSSAVG).  The power
%     spectra is plotted on a map with a Robinson projection and the map
%     limits are scaled to fit the data & station positions.  Note that
%     this function requires the data positions to be sampled in a regular
%     grid (see MESHGRID).  The default plots the coastlines at a crude
%     resolution.  Additional options to change the map are allowed through
%     the final Usage form.
%
%     PLOTGEOFSS(S,DBLIM) sets the dB limits for the power spectra
%     colormap.  The default is [-12 0] for the default ZERODB (see the
%     next Usage form).  If ZERODB IS 'min' or 'median', the default DBLIM
%     is [0 12].  DBLIM must be a real-valued 2-element vector.
%
%     PLOTGEOFSS(S,DBLIM,ZERODB) changes what 0dB corresponds to in the
%     plot.  The allowed values are 'min', 'max', 'median', & 'abs'.  The
%     default is 'max'.
%
%     PLOTGEOFSS(S,DBLIM,ZERODB,'PARAM1',VAL1,...) passes additional
%     parameter values for mapping.  See MMAP for details.
%
%    Notes:
%
%    Examples:
%     % In search of the 26s microseism:
%     [lat,lon]=meshgrid(-60:60,-60:60);
%     s=geofssxc(xcdata,[lat(:) lon(:)],27:33,[1/26.3 1/26]);
%     s0=geofssavg(s);
%     plotgeofss(s0);
%
%    See also: GEOFSS, GEOFSSXC, GEOFSSHORZXC, GEOFSSARF, GEOFSSAVG,
%              GEOFSSSUB, GEOFSSINFO, GEOFSSCORR, GEOFSSFRAMESLIDE,
%              GEOFSSFREQSLIDE, GEOFSSSLOWSLIDE, PLOTGEOFSSARF

%     Version History:
%        June 25, 2010 - initial version
%        July  1, 2010 - no land cover
%        July  6, 2010 - update for new struct
%        July 14, 2010 - added some lines for not plotting stations
%        Oct. 10, 2010 - all plotting functions use proper ax calls, tagged
%                        plots as 'fkmap'
%        Dec.  8, 2010 - use '^o' for deg symbol rather than \circ
%        Feb. 16, 2011 - color code fix
%        Apr.  4, 2012 - minor doc update
%        Apr. 25, 2012 - use nanmedian for median determination
%        June  4, 2012 - altered from plotgeofkmap, full mmap options
%        June 10, 2012 - handle full method
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 10, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));
if(nargin>3 && ~mod(nargin,2))
    error('seizmo:plotgeofss:badInput',...
        'Unpaired MMAP parameter/value!');
end

% check geofss struct
error(chkgeofss(s));

% require scalar struct
if(~isscalar(s))
    error('seizmo:plotgeofss:badInput',...
        'S must be a scalar geofss struct');
end

% require scalar freq & slow dimensions
if(any(s.vector))
    error('seizmo:plotgeofss:badInput',...
        'S needs to be reduced using GEOFSSAVG!');
end

% default/check scaling type
if(nargin<3 || isempty(zerodb)); zerodb='max'; end
if(~ischar(zerodb) ...
        || ~ismember(lower(zerodb),{'max' 'min' 'median' 'med' 'abs'}))
    error('seizmo:plotgeofss:badSTYPE',...
        'STYPE must be ''min'' ''max'' ''median'' or ''abs''!');
end
zerodb=lower(zerodb);

% default/check dblim
if(nargin<2 || isempty(dblim))
    switch zerodb
        case {'min' 'median' 'med'}
            dblim=[0 12];
        case {'max' 'abs'}
            dblim=[-12 0];
    end
end
if(~isreal(dblim) || numel(dblim)~=2)
    error('seizmo:plotgeofss:badDBLIM',...
        'DBLIM must be a real valued 2 element vector!');
end
dblim=sort([dblim(1) dblim(2)]);

% convert to dB
switch s.method
    case {'coarray' 'full'}
        s.spectra=10*log10(abs(real(s.spectra)));
    otherwise
        s.spectra=10*log10(s.spectra);
end

% rescale response
switch zerodb
    case 'min'
        zdb=min(s.spectra(:));
        s.spectra=s.spectra-zdb;
    case 'max'
        zdb=max(s.spectra(:));
        s.spectra=s.spectra-zdb;
    case {'median' 'med'}
        zdb=nanmedian(s.spectra(:));
        s.spectra=s.spectra-zdb;
    case 'abs'
        zdb=0;
end

% reshape spectra for plotting & to account for pcolor
nlat=numel(unique(s.latlon(:,1)));
nlon=numel(unique(s.latlon(:,2)));
s.latlon=reshape(s.latlon,[nlon nlat 2]); % this breaks if non-meshgrid
s.spectra=reshape(s.spectra,[nlon nlat]); % this breaks if non-meshgrid
latstep=s.latlon(1,2,1)-s.latlon(1,1,1);
lonstep=s.latlon(2,1,2)-s.latlon(1,1,2);
s.latlon(:,:,1)=s.latlon(:,:,1)-latstep/2;
s.latlon(:,:,2)=s.latlon(:,:,2)-lonstep/2;
s.latlon=s.latlon([1:end end],[1:end end],:);
s.latlon(:,end,1)=s.latlon(:,end,1)+latstep;
s.latlon(end,:,2)=s.latlon(end,:,2)+lonstep;
s.spectra=s.spectra([1:end end],[1:end end]);

% limits of map
minlat=min([s.stla; min(s.latlon(:,:,1))']);
maxlat=max([s.stla; max(s.latlon(:,:,1))']);
minlon=min([s.stlo; min(s.latlon(:,:,2))']);
maxlon=max([s.stlo; max(s.latlon(:,:,2))']);

% plot map
ax=mmap('po',{'lat',[minlat maxlat],'lon',[minlon maxlon]},...
    'land',false,'sea',false,varargin{:});

% access to m_map globals for map boundaries
global MAP_VAR_LIST

% plot spectra (up to 3 times)
hold(ax,'on');
if(any(s.latlon(:,:,2)>=MAP_VAR_LIST.longs(1) ...
        & s.latlon(:,:,2)<=MAP_VAR_LIST.longs(2)))
    m_pcolor(s.latlon(:,:,2),s.latlon(:,:,1),double(s.spectra),...
        'parent',ax);
end
if(any(s.latlon(:,:,2)-360>=MAP_VAR_LIST.longs(1) ...
        & s.latlon(:,:,2)-360<=MAP_VAR_LIST.longs(2)))
    m_pcolor(s.latlon(:,:,2)-360,s.latlon(:,:,1),double(s.spectra),...
        'parent',ax);
end
if(any(s.latlon(:,:,2)+360>=MAP_VAR_LIST.longs(1) ...
        & s.latlon(:,:,2)+360<=MAP_VAR_LIST.longs(2)))
    m_pcolor(s.latlon(:,:,2)+360,s.latlon(:,:,1),double(s.spectra),...
        'parent',ax);
end

% move spectra back
movekids(findobj(ax,'tag','m_pcolor'),'back');
movekids(findobj(ax,'tag','m_grid_color'),'back');

% this hides some warnings
ocean=get(findobj(ax,'tag','m_grid_color'),'facecolor');
if(~ischar(ocean))
    set(findobj(ax,'tag','m_grid_color'),'facevertexcdata',ocean(ones(...
        size(get(findobj(ax,'tag','m_grid_color'),'vertices'),1),1),:));
else
    delete(findobj(ax,'tag','m_grid_color'));
end

% modify
shading(ax,'flat');
bg=get(get(ax,'parent'),'color');
fg=get(findobj(ax,'tag','m_grid_box'),'color');
if(strcmp(bg,'w') || isequal(bg,[1 1 1]))
    colormap(ax,flipud(fire));
elseif(strcmp(bg,'k') || isequal(bg,[0 0 0]))
    colormap(ax,fire);
end
set(ax,'clim',dblim);
hold(ax,'off');

% wrap station longitudes to within 180deg of plot center
while(any(abs(s.stlo-mean(MAP_VAR_LIST.longs))>180))
    s.stlo(s.stlo<MAP_VAR_LIST.longs(1))=...
        s.stlo(s.stlo<MAP_VAR_LIST.longs(1))+360;
    s.stlo(s.stlo>MAP_VAR_LIST.longs(2))=...
        s.stlo(s.stlo>MAP_VAR_LIST.longs(2))-360;
end

% map the stations
% - have to do this afterwards b/c of matlab warnings/issues
hold(ax,'on');
mmap('st',[s.stla s.stlo],'parent',ax,varargin{:});
hold(ax,'off');

% colorbar & title
c=colorbar('eastoutside','peer',ax,'xcolor',fg,'ycolor',fg);
xlabel(c,'dB','color',fg);
fmin=min(s.freq); fmax=max(s.freq);
smin=min(s.slow); smax=max(s.slow);
title(ax,{[] ['Number of Stations:  ' num2str(s.nsta)] ...
    ['Begin Time:  ' sprintf('%d.%03d %02d:%02d:%02g',s.butc) ' UTC'] ...
    ['End Time  :  ' sprintf('%d.%03d %02d:%02d:%02g',s.eutc) ' UTC'] ...
    ['Period    :  ' num2str(1/fmax) ' to ' num2str(1/fmin) ' s'] ...
    ['Slowness  :  ' num2str(smin) ' to ' num2str(smax) ' s/^o'] ...
    ['0 dB = ' num2str(zdb) 'dB'] []},'color',fg);

% set zerodb in userdata
% - this is for plotgeofssupdate
userdata.zerodb=zerodb;
set(ax,'userdata',userdata,'createfcn','plotgeofss');

% return figure handle
if(nargout); varargout{1}=ax; end

end
