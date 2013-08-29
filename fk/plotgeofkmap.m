function [varargout]=plotgeofkmap(map,dblim,zerodb,varargin)
%PLOTGEOFKMAP    Plots frequency-slowness-position beam data
%
%    Usage:    plotgeofkmap(map)
%              plotgeofkmap(map,dblim)
%              plotgeofkmap(map,dblim,zerodb)
%              plotgeofkmap(map,dblim,zerodb,'mmap_opt1',mmap_val1,...)
%              ax=plotgeofkmap(...)
%
%    Description:
%     PLOTGEOFKMAP(MAP) plots the frequency-slowness-position beam data in
%     geofk struct MAP.  See a geofk function like GEOFKXCVOLUME for
%     details on the struct.  The data is plotted on a map with a
%     Robinson projection and the map limits are scaled to fit the
%     beam data & station positions.  Note that the beam data positions
%     should form a regular grid (using a function like MESHGRID).
%
%     PLOTGEOFKMAP(MAP,DBLIM) sets the dB limits for coloring the
%     response info.  The default is [-12 0] for the default ZERODB (see
%     next Usage form).  If ZERODB IS 'min' or 'median', the default DBLIM
%     is [0 12].  DBLIM must be a real-valued 2-element vector.
%
%     PLOTGEOFKMAP(MAP,DBLIM,ZERODB) changes what 0dB corresponds
%     to in the plot.  The allowed values are 'min', 'max', 'median', &
%     'abs'.  The default is 'max'.
%
%     PLOTGEOFKMAP(MAP,DBLIM,ZERODB,'MMAP_OPT1',MMAP_VAL1,...) passes
%     additional options on to MMAP to alter the map.
%
%     AX=PLOTGEOFKMAP(...) returns the axes drawn in.
%
%    Notes:
%
%    Examples:
%     % In search of the 26s microseism:
%     [lat,lon]=meshgrid(-60:60,-60:60);
%     zgeo=geofkxcvolume(xcdata,[lat(:) lon(:)],27:33,[1/26.3 1/26]);
%     zgeo0=geofkvol2map(zgeo);
%     plotgeofkmap(zgeo0);
%
%    See also: GEOFKFREQSLIDE, GEOFKSLOWSLIDE, GEOFKVOL2MAP, GEOFKXCVOLUME,
%              GEOFKXCHORZVOLUME, CHKGEOFKSTRUCT, UPDATEGEOFKMAP

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
%        Aug. 15, 2012 - plot coasts and latlon grid only
%        Aug. 28, 2013 - use mmap image option, reorder inputs
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 28, 2013 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% check fk struct
error(chkgeofkstruct(map));

% don't allow array/volume
if(~isscalar(map) || any(map.volume))
    error('seizmo:plotgeofkmap:badInput',...
        'MAP must be a scalar geofk struct and not a volume!');
end

% default/check scaling type
if(nargin<3 || isempty(zerodb)); zerodb='max'; end
if(~ischar(zerodb) ...
        || ~ismember(lower(zerodb),{'max' 'min' 'median' 'abs'}))
    error('seizmo:plotgeofkmap:badSTYPE',...
        'STYPE must be ''min'' ''max'' ''median'' or ''abs''!');
end
zerodb=lower(zerodb);

% default/check dblim
if(nargin<2 || isempty(dblim))
    switch zerodb
        case {'min' 'median'}
            dblim=[0 12];
        case {'max' 'abs'}
            dblim=[-12 0];
    end
end
if(~isreal(dblim) || numel(dblim)~=2)
    error('seizmo:plotgeofkmap:badDBLIM',...
        'DBLIM must be a real valued 2 element vector!');
end
dblim=sort([dblim(1) dblim(2)]);

% rescale response
switch zerodb
    case 'min'
        map.normdb=map.normdb+min(map.beam(:));
        map.beam=map.beam-min(map.beam(:));
    case 'max'
        map.normdb=map.normdb+max(map.beam(:));
        map.beam=map.beam-max(map.beam(:));
    case 'median'
        map.normdb=map.normdb+nanmedian(map.beam(:));
        map.beam=map.beam-nanmedian(map.beam(:));
    case 'abs'
        map.beam=map.beam+map.normdb;
        map.normdb=0;
end

% reshape beam & account for pcolor
nlat=numel(unique(map.latlon(:,1)));
nlon=numel(unique(map.latlon(:,2)));
map.latlon=reshape(map.latlon,[nlon nlat 2]);
map.beam=reshape(map.beam,[nlon nlat]);

% get max/min lat/lon of map & stations
minlat=min([map.stla; min(map.latlon(:,:,1))']);
maxlat=max([map.stla; max(map.latlon(:,:,1))']);
minlon=min([map.stlo; min(map.latlon(:,:,2))']);
maxlon=max([map.stlo; max(map.latlon(:,:,2))']);

% a couple mmap default changes
% - use min/max of lat/lon as the map boundary
% - do not show land/ocean
varargin=[{'po' {'lat' [minlat maxlat] ...
    'lon' [minlon maxlon]} 'l' false 'o' false} varargin];

% draw map
ax=mmap('image',{map.latlon(:,:,1) map.latlon(:,:,2) double(map.beam)},...
    'st',[map.stla map.stlo],varargin{:});

% extract color
bg=get(get(ax,'parent'),'color');
fg=get(findobj(ax,'tag','m_grid_box'),'color');

% modify
if(strcmp(bg,'w') || isequal(bg,[1 1 1]))
    colormap(ax,flipud(fire));
elseif(strcmp(bg,'k') || isequal(bg,[0 0 0]))
    colormap(ax,fire);
else
    if(ischar(bg)); bg=name2rgb(bg); end
    hsv=rgb2hsv(bg);
    colormap(ax,hsvcustom(hsv));
end
set(ax,'clim',dblim);

% colorbar & title
c=colorbar('eastoutside','peer',ax,'xcolor',fg,'ycolor',fg);
xlabel(c,'dB','color',fg);
fmin=min(map.freq); fmax=max(map.freq);
smn=min(map.horzslow); smx=max(map.horzslow);
title(ax,{[] ['Number of Stations:  ' num2str(map.nsta)] ...
    ['Begin Time:  ' sprintf('%d.%03d %02d:%02d:%02g',map.butc) ' UTC'] ...
    ['End Time  :  ' sprintf('%d.%03d %02d:%02d:%02g',map.eutc) ' UTC'] ...
    ['Period    :  ' num2str(1/fmax) ' to ' num2str(1/fmin) ' s'] ...
    ['Horiz. Slowness :  ' num2str(smn) ' to ' num2str(smx) ' s/^o'] ...
    ['0 dB = ' num2str(map.normdb) 'dB'] []},'color',fg);

% set zerodb & dblim in userdata
% - this is for updategeofkmap
userdata.zerodb=zerodb;
userdata.dblim=dblim;
set(ax,'userdata',userdata,'tag','geofkmap');

% return figure handle
if(nargout); varargout{1}=ax; end

end
