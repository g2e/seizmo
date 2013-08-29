function [varargout]=plotbathyexcite(c,la,lo,crng,varargin)
%PLOTBATHYEXCITE    Plots 2ndary microseism bathymetric excitation coeff
%
%    Usage:    plotbathyexcite(c,lat,lon)
%              plotbathyexcite(c,lat,lon,clim)
%              plotbathyexcite(c,lat,lon,clim,'mmap_opt1',mmap_val1,...)
%              ax=plotbathyexcite(...)
%
%    Description:
%     PLOTBATHYEXCITE(C,LAT,LON) maps bathymetric coefficient data given
%     in the 2D array C.  The data should be regularly sampled in latitude
%     and longitude (given by vectors LAT & LON).  Dimensions of C should
%     be NLATxNLON.
%
%     PLOTBATHYEXCITE(C,LAT,LON,CLIM) sets the colormap limits of the
%     coefficients.  The default is [0 1].  CLIM must be a real-valued
%     2-element vector.
%     
%     PLOTBATHYEXCITE(C,LAT,LON,CLIM,'MMAP_OPT1',MMAP_VAL1,...) passes
%     additional options on to MMAP to alter the map.
%
%     AX=PLOTBATHYEXCITE(...) returns the axes drawn in.
%
%    Notes:
%
%    Examples:
%     % Get bathymetric excitation coefficients for Crust2.0 and plot:
%     [lon,lat]=meshgrid(-179:2:179,89:-2:-89);
%     c2elev=getc2elev(lat,lon);
%     c2elev(c2elev>0)=0; % mask out land
%     c2=getcrust2(lat,lon);
%     vs=cat(1,c2.vs);
%     vs=reshape(vs(:,3:7),90,180,5);
%     thick=cat(1,c2.thick);
%     thick=reshape(thick(:,3:7),90,180,5);
%     c2vsavg=sum(thick,3)./sum(thick./vs,3);
%     c=bathy_micro_excite(-c2elev,1/7.5,c2vsavg*1000);
%     ax=plotbathyexcite(c,lat,lon);
%     title(ax,{[] 'Crust2.0 Bathymetric Excitation Coefficient Map' ...
%         'Period: 7.5s   Vs: 2.2-3.75km/s' []},'color','w');
%
%    See also: BATHY_MICRO_EXCITE

%     Version History:
%        Feb. 15, 2011 - initial version
%        May   5, 2012 - minor doc update
%        May  18, 2012 - improved label of the colorbar, use gray colormap
%        Aug. 27, 2013 - use mmap image option
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 27, 2013 at 15:05 GMT

% todo:

% check nargin
error(nargchk(3,inf,nargin));

% check coefficient array
if(~isreal(c) || any(c(:)<0) || ndims(c)~=2)
    error('seizmo:plotbathyexcite:badInput',...
        'C must be a 2D array of real positive values!');
end

% check lat/lon match up
if(~isreal(la) || ~isreal(lo) || any(abs(la(:))>90))
    error('seizmo:plotbathyexcite:badInput',...
        'LAT & LON must be real-valued arrays in the appropriate range!');
end
if(~isvector(la)); la=flipud(unique(la(:))); end
if(~isvector(lo)); lo=unique(lo(:)); end
if(~isequal(size(c),[numel(la) numel(lo)]))
    error('seizmo:plotbathyexcite:badInput',...
        'LAT & LON must be vectors matching C dimensions!');
end

% default/check crng
if(nargin<4 || isempty(crng)); crng=[0 1]; end
if(~isreal(crng) || numel(crng)~=2)
    error('seizmo:plotbathyexcite:badInput',...
        'CRNG must be a real valued 2 element vector!');
end
crng=sort([crng(1) crng(2)]);

% a couple mmap default changes
% - use min/max of lat/lon as the map boundary
% - do not show land/ocean
varargin=[{'po' {'lat' [min(la) max(la)] ...
    'lon' [min(lo) max(lo)]} 'l' false 'o' false} varargin];

% draw map
ax=mmap('image',{la lo c.'},varargin{:});

% extract color
bg=get(get(ax,'parent'),'color');
fg=get(findobj(ax,'tag','m_grid_box'),'color');

% set colormap
if(strcmp(bg,'w') || isequal(bg,[1 1 1]))
    colormap(ax,flipud(gray));
elseif(strcmp(bg,'k') || isequal(bg,[0 0 0]))
    colormap(ax,gray);
else
    if(ischar(bg)); bg=name2rgb(bg); end
    hsv=rgb2hsv(bg);
    colormap(ax,hsvcustom(hsv));
end
set(ax,'clim',crng);

% colorbar & title
c=colorbar('eastoutside','peer',ax,'xcolor',fg,'ycolor',fg);
xlabel(c,'$$\sum {c_{m}^2}$$','color',fg,'interpreter','latex');
title(ax,{[] 'Bathymetric Excitation Coefficient Map' []},'color',fg);

% return figure handle
set(ax,'tag','bathyexcitemap');
if(nargout); varargout{1}=ax; end

end
