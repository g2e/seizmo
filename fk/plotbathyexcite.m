function [varargout]=plotbathyexcite(c,la,lo,crng,popt,fgcolor,bgcolor,ax)
%PLOTBATHYEXCITE    Plots 2ndary microseism bathymetric excitation coeff
%
%    Usage:    plotbathyexcite(c,lat,lon)
%              plotbathyexcite(c,lat,lon,clim)
%              plotbathyexcite(c,lat,lon,clim,projopt)
%              plotbathyexcite(c,lat,lon,clim,projopt,fgcolor,bgcolor)
%              plotbathyexcite(c,lat,lon,clim,projopt,fgcolor,bgcolor,ax)
%              ax=plotbathyexcite(...)
%
%    Description:
%     PLOTBATHYEXCITE(C,LAT,LON) plots bathymetric coefficient data given
%     in the 2D array C.  The data should be regularly sampled in latitude
%     and longitude (given by vectors LAT & LON).  The data is plotted on a
%     map with a Hammer-Aitoff projection and the map limits are scaled to
%     fit the latitude and longitude limits.  This plots GSHHS coastlines
%     and borders in crude-resolution which may take a few moments - please
%     be patient. Dimensions of C should be NUMEL(LAT)xNUMEL(LON).
%
%     PLOTBATHYEXCITE(C,LAT,LON,CLIM) sets the coloring limits of the
%     coefficients limits for coloring.  The default is [0 1].  CLIM must
%     be a real-valued 2-element vector.
%
%     PLOTBATHYEXCITE(C,LAT,LON,CLIM,PROJOPT) allows passing options to
%     M_PROJ.  See M_PROJ('SET') for possible projections and see
%     M_PROJ('GET',PROJ) for a list of possible additional options specific
%     to that projection.
%
%     PLOTBATHYEXCITE(C,LAT,LON,CLIM,PROJOPT,FGCOLOR,BGCOLOR) specifies
%     foreground and background colors of the plot.  The default is 'w' for
%     FGCOLOR & 'k' for BGCOLOR.  Note that if one is specified and the
%     other is not, an opposing color is found using INVERTCOLOR.  The
%     color scale is also changed so the lower color clip is at BGCOLOR.
%
%     PLOTBATHYEXCITE(MAP,PROJOPT,DBLIM,ZERODB,FGCOLOR,BGCOLOR,AX) sets the
%     axes to draw in.  This is useful for subplots, guis, etc.
%
%     AX=PLOTBATHYEXCITE(...)
%
%    Notes:
%     - This function can be easily altered to work for any image data.
%
%    Examples:
%     % Get bathymetric excitation coefficients for Crust2.0 and plot:
%     [lon,lat]=meshgrid(-179:2:179,89:-2:-89);
%     c2elev=getc2elev(lat,lon);
%     c2elev(c2elev>0)=0; % mask out land
%     c=bathy_micro_excite(-c2elev);
%     ax=plotbathyexcite(c,lat,lon);
%     title(ax,{[] 'Crust2.0 Bathymetric Excitation Coefficient Map' ...
%         'Period: 6.28s   Vs: 2.8km/s' []},'color','w');
%
%    See also: BATHY_MICRO_EXCITE

%     Version History:
%        Feb. 15, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 15, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,8,nargin));

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

% check colors
if(nargin<6);
    fgcolor='w'; bgcolor='k';
elseif(nargin<7)
    if(isempty(fgcolor))
        fgcolor='w'; bgcolor='k';
    else
        bgcolor=invertcolor(fgcolor,true);
    end
else
    if(isempty(fgcolor))
        if(isempty(bgcolor))
            fgcolor='w'; bgcolor='k';
        else
            fgcolor=invertcolor(bgcolor,true);
        end
    elseif(isempty(bgcolor))
        bgcolor=invertcolor(fgcolor,true);
    end
end

% change char to something rgb
if(ischar(fgcolor)); fgcolor=name2rgb(fgcolor); end
if(ischar(bgcolor)); bgcolor=name2rgb(bgcolor); end

% check handle
if(nargin<8 || isempty(ax) || ~isscalar(ax) || ~isreal(ax) ...
        || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
    figure('color',bgcolor);
    ax=gca;
else
    axes(ax);
    h=get(ax,'children'); delete(h);
    h=findobj(get(get(ax,'parent'),'children'),'peer',ax); delete(h);
end

% map colors & coast/border res
gshhs='c';
ocean=bgcolor;
land=fgcolor;
border=fgcolor;

% access to m_map globals for map boundaries
global MAP_VAR_LIST

% reshape beam & account for pcolor
latstep=la(2)-la(1);
lonstep=lo(2)-lo(1);
la=[la(:)-latstep/2; la(end)+latstep/2];
lo=[lo(:)'-lonstep/2 lo(end)+lonstep/2];
c=c([1:end end],[1:end end]);

% get max/min lat/lon of map & stations
minlat=min(la); maxlat=max(la);
minlon=min(lo); maxlon=max(lo);

% default/check projopt
if(nargin<5 || isempty(popt))
    popt={'robinson','lat',[minlat maxlat],'lon',[minlon maxlon]};
end
if(ischar(popt)); popt=cellstr(popt); end
if(~iscell(popt))
    error('seizmo:plotbathyexcite:badInput',...
        'PROJOPT must be a cell array of args for M_PROJ!');
end

% setup projection
m_proj(popt{:});
set(ax,'color',ocean);

% plot geofk beam
hold(ax,'on');
if(any(lo>MAP_VAR_LIST.longs(1) & lo<MAP_VAR_LIST.longs(2)))
    m_pcolor(lo,la,c,'parent',ax);
end
if(any(lo-360>MAP_VAR_LIST.longs(1) & lo-360<MAP_VAR_LIST.longs(2)))
    m_pcolor(lo-360,la,c,'parent',ax);
end
if(any(lo+360>MAP_VAR_LIST.longs(1) & lo+360<MAP_VAR_LIST.longs(2)))
    m_pcolor(lo+360,la,c,'parent',ax);
end

% modify
shading(ax,'flat');
if(strcmp(bgcolor,'w') || isequal(bgcolor,[1 1 1]))
    colormap(ax,flipud(fire));
elseif(strcmp(bgcolor,'k') || isequal(bgcolor,[0 0 0]))
    colormap(ax,fire);
else
    if(ischar(bgcolor))
        bgcolor=name2rgb(bgcolor);
    end
    hsv=rgb2hsv(bgcolor);
    colormap(ax,hsvcustom(hsv));
end
set(ax,'clim',crng);
hold(ax,'off');

% now add coastlines and political boundaries
axes(ax);
m_gshhs([gshhs 'c'],'color',land);
m_gshhs([gshhs 'b'],'color',border);
m_grid('color',fgcolor);

% hackery to color oceans at large when the above fails
set(findobj(ax,'tag','m_grid_color'),'facecolor',ocean);

% colorbar & title
c=colorbar('eastoutside','peer',ax,'xcolor',fgcolor,'ycolor',fgcolor);
xlabel(c,'c','color',fgcolor);
title(ax,{[] 'Bathymetric Excitation Coefficient Map' []},...
    'color',fgcolor);

% tag plot
set(ax,'tag','bathyexcite');

% return figure handle
if(nargout); varargout{1}=ax; end

end
