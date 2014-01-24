function [varargout]=plotbathyexcite(ph,c,la,lo,clim,varargin)
%PLOTBATHYEXCITE    Plots 2ndary microseism bathymetric excitation coeff
%
%    Usage:    plotbathyexcite(ph,c,lat,lon)
%              plotbathyexcite(ph,c,lat,lon,clim)
%              plotbathyexcite(ph,c,lat,lon,clim,'mmap_opt1',mmap_val1,...)
%              ax=plotbathyexcite(...)
%
%    Description:
%     PLOTBATHYEXCITE(PH,C,LAT,LON) maps bathymetric coefficient data given
%     in the 2D array C.  The data should be regularly sampled in latitude
%     and longitude (given by vectors LAT & LON).  Dimensions of C should
%     be NLATxNLON.  The output map will display the squared values of C.
%     PH is the phase string (such as 'P') and is used for labeling.
%
%     PLOTBATHYEXCITE(PH,C,LAT,LON,CLIM) sets the colormap limits of the
%     coefficients.  The default is [0 1] which is fine for Rayleigh wave
%     excitation coefficients but not for P or S waves.  CLIM must be a
%     real-valued 2-element vector.
%     
%     PLOTBATHYEXCITE(PH,C,LAT,LON,CLIM,'MMAP_OPT1',MMAP_VAL1,...) passes
%     additional options on to MMAP to alter the map.
%
%     AX=PLOTBATHYEXCITE(...) returns the axes drawn in.
%
%    Notes:
%
%    Examples:
%     % Rayleigh Bathymetric Excitation Map:
%     [lon,lat]=meshgrid(-179.5:179.5,-89.5:89.5);
%     mod=getcrust(lat,lon);
%     mod.top(mod.top(:,2)>0,2)=0; % mask out land
%     mod.vs(mod.vs==0)=1; % avoid divide by zero
%     vsavg=sum(mod.thk(:,2:8),2)./sum(mod.thk(:,2:8)./mod.vs(:,2:8),2);
%     c=bathy_micro_excite('R',-mod.top(:,2)*1000,1/7.5,vsavg*1000);
%     ax=plotbathyexcite('R',reshape(c,[180 360]),lat,lon);
%     title(ax,{[] ...
%         'Rayleigh Bathymetric Excitation Map for' ...
%         'Wave-Wave Interference of Ocean Waves' ...
%         'Seismic Period: 7.5s' []},'color','w');
%
%    See also: BATHY_MICRO_EXCITE

%     Version History:
%        Feb. 15, 2011 - initial version
%        May   5, 2012 - minor doc update
%        May  18, 2012 - improved label of the colorbar, use gray colormap
%        Aug. 27, 2013 - use mmap image option
%        Jan. 14, 2014 - bugfix: squares input c values, added phase string
%                        argument (required), changed crng to clim
%        Jan. 23, 2014 - doc update, made phase argument first to match
%                        bathy_micro_excite, fix clim bug, no title set
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 23, 2014 at 15:05 GMT

% todo:

% check nargin
error(nargchk(4,inf,nargin));

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
if(~isvector(la))
    if(diff(la(1:2)))
        % lats down the column
        la=la(:,1);
    else
        % lats across the row
        la=la(1,:);
    end
end
if(~isvector(lo))
    if(diff(lo(1:2)))
        % lons down the column
        lo=lo(:,1);
    else
        % lons across the row
        lo=lo(1,:);
    end
if(~isequal(size(c),[numel(la) numel(lo)]))
    error('seizmo:plotbathyexcite:badInput',...
        'LAT & LON must be vectors matching C dimensions!');
end

% require ph is a string
if(~ischar(ph) || ndims(ph)>2 || size(ph,1)>1)
    error('seizmo:plotbathyexcite:badInput',...
        'PH must be a string!');
end

% default/check clim
if(nargin<5 || isempty(clim)); clim=[0 1]; end
if(~isreal(clim) || numel(clim)~=2)
    error('seizmo:plotbathyexcite:badInput',...
        'CLIM must be a real valued 2 element vector!');
end
clim=sort([clim(1) clim(2)]);

% a couple mmap default changes
% - use min/max of lat/lon as the map boundary
% - do not show land/ocean
varargin=[{'po' {'lat' [min(la) max(la)] ...
    'lon' [min(lo) max(lo)]} 'l' false 'o' false} varargin];

% square input coefficients to get proper amplification
c=c.^2;

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
set(ax,'clim',clim);

% colorbar
c=colorbar('eastoutside','peer',ax,'xcolor',fg,'ycolor',fg);
ph=['$$c_{' ph '}^2$$'];
xlabel(c,ph,'color',fg,'interpreter','latex');

% return figure handle
set(ax,'tag','bathyexcitemap');
if(nargout); varargout{1}=ax; end

end
