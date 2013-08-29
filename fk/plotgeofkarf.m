function [varargout]=plotgeofkarf(arf,dblim,zerodb,varargin)
%PLOTGEOFKARF    Plots geofk array response
%
%    Usage:    plotgeofkarf(arf)
%              plotgeofkarf(arf,dblim)
%              plotgeofkarf(arf,dblim,zerodb)
%              plotgeofkarf(arf,dblim,zerodb,'mmap_opt1',mmap_val1,...)
%              ax=plotgeofkarf(...)
%
%    Description:
%     PLOTGEOFKARF(ARF) plots the frequency-slowness-position beam data in
%     geofk struct ARF.  See a geofk function like GEOFKARF for
%     details on the struct.  The data is plotted on a map with a
%     Robinson projection and the map limits are scaled to fit the
%     beam data & station positions.  Note that the beam data positions
%     should form a regular grid (using a function like MESHGRID).
%
%     PLOTGEOFKARF(ARF,DBLIM) sets the dB limits for coloring the
%     response info.  The default is [-12 0] for the default ZERODB (see
%     next Usage form).  If ZERODB IS 'min' or 'median', the default DBLIM
%     is [0 12].  DBLIM must be a real-valued 2-element vector.
%
%     PLOTGEOFKARF(ARF,DBLIM,ZERODB) changes what 0dB corresponds
%     to in the plot.  The allowed values are 'min', 'max', 'median', &
%     'abs'.  The default is 'max'.
%
%     PLOTGEOFKARF(ARF,DBLIM,ZERODB,'MMAP_OPT1',MMAP_VAL1,...) passes
%     additional options on to MMAP to alter the map.
%
%     AX=PLOTGEOFKARF(...) returns the axes drawn in.
%
%    Notes:
%
%    Examples:
%     % ARF at 500s for a 30 station global array for 6 sources:
%     st=randlatlon(30);
%     ev=randlatlon(6);
%     [la,lo]=meshgrid(-90:90,-180:180);
%     arf=geofkarf(st,[la(:) lo(:)],30,ev,30,1/500,'center');
%     arf=geofkarf2map(arf);
%     plotgeofkarf(arf,[-6 0]);
%
%    See also: GEOFKARF, GEOFKARF2MAP, GEOFKSUBARF, UPDATEGEOFKARF,
%              GEOFKARFSLOWSLIDE, CHKGEOFKARFSTRUCT

%     Version History:
%        July  7, 2010 - update for new struct
%        Oct.  6, 2010 - truncate title if too many ARF locations
%        Oct. 10, 2010 - all plotting functions use proper ax calls, tagged
%                        plots as 'fkmap'
%        Dec.  8, 2010 - use '^o' for deg symbol rather than \circ
%        Feb. 16, 2011 - color code fix
%        Feb.  2, 2012 - use robinson projection like plotgeofkmap
%        Apr.  4, 2012 - minor doc update
%        May   5, 2012 - minor doc update
%        Aug. 28, 2013 - use mmap image option, reorder inputs, add example
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 28, 2013 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% check fk struct
error(chkgeofkarfstruct(arf));

% don't allow array/volume
if(~isscalar(arf) || any(arf.volume(2)))
    error('seizmo:plotgeofkarf:badInput',...
        'ARF must be a scalar geofk struct and not a volume!');
end

% default/check scaling type
if(nargin<3 || isempty(zerodb)); zerodb='max'; end
if(~ischar(zerodb) ...
        || ~ismember(lower(zerodb),{'max' 'min' 'median' 'abs'}))
    error('seizmo:plotgeofkarf:badSTYPE',...
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
    error('seizmo:plotgeofkarf:badDBLIM',...
        'DBLIM must be a real valued 2 element vector!');
end
dblim=sort([dblim(1) dblim(2)]);

% rescale response
switch zerodb
    case 'min'
        arf.normdb=arf.normdb+min(arf.beam(:));
        arf.beam=arf.beam-min(arf.beam(:));
    case 'max'
        arf.normdb=arf.normdb+max(arf.beam(:));
        arf.beam=arf.beam-max(arf.beam(:));
    case 'median'
        arf.normdb=arf.normdb+median(arf.beam(:));
        arf.beam=arf.beam-median(arf.beam(:));
    case 'abs'
        arf.beam=arf.beam+arf.normdb;
        arf.normdb=0;
end

% reshape beam
nlat=numel(unique(arf.latlon(:,1)));
nlon=numel(unique(arf.latlon(:,2)));
arf.latlon=reshape(arf.latlon,[nlon nlat 2]);
arf.beam=reshape(arf.beam,[nlon nlat]);

% get max/min lat/lon of arf & stations
minlat=min([arf.stla; min(arf.latlon(:,:,1))']);
maxlat=max([arf.stla; max(arf.latlon(:,:,1))']);
minlon=min([arf.stlo; min(arf.latlon(:,:,2))']);
maxlon=max([arf.stlo; max(arf.latlon(:,:,2))']);

% a couple mmap default changes
% - use min/max of lat/lon as the map boundary
% - do not show land/ocean
varargin=[{'po' {'lat' [minlat maxlat] ...
    'lon' [minlon maxlon]} 'l' false 'o' false} varargin];

% draw map
ax=mmap('image',{arf.latlon(:,:,1) arf.latlon(:,:,2) double(arf.beam)},...
    'st',[arf.stla arf.stlo],'ev',arf.latlon0,varargin{:});

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
smn=min(arf.horzslow); smx=max(arf.horzslow);
if(arf.nsw<=5)
    titstr=cell(arf.nsw,1);
    for i=1:arf.nsw
        titstr{i}=sprintf(['SLOWNESS: %gs/^o, LAT: %g^o, ' ...
            'LON: %g^o, PERIOD: %gs'],arf.horzslow0(i),...
            arf.latlon0(i,1),arf.latlon0(i,2),1/arf.freq0(i));
    end
else
    titstr{1}=[num2str(arf.nsw) ' Locations'];
end
title(ax,[{[]}; 'Array Response Function @ '; titstr; ...
    ['Number of Stations: ' num2str(arf.nsta)]; ...
    ['Horiz. Slowness : ' num2str(smn) ' to ' num2str(smx) ' s/^o']; ...
    ['0 dB = ' num2str(arf.normdb) 'dB']; {[]}],'color',fg);

% set zerodb & dblim in userdata
% - this is for updategeofkarf
userdata.zerodb=zerodb;
userdata.dblim=dblim;
set(ax,'userdata',userdata,'tag','geofkarf');

% return figure handle
if(nargout); varargout{1}=ax; end

end
