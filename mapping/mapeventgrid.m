function []=mapeventgrid(ax,evla,evlo,c,rrmajor,rrminor,azmajor,azminor)
%MAPEVENTGRID    Draws range/azimuth grid (relative to an event) on a map
%
%    Usage:    mapeventgrid(ax,evla,evlo)
%              mapeventgrid(ax,evla,evlo,c)
%              mapeventgrid(ax,evla,evlo,c,rrmaj,rrmin)
%              mapeventgrid(ax,evla,evlo,c,rrmaj,rrmin,azmaj,azmin)
%
%    Description:
%     MAPEVENTGRID(AX,EVLA,EVLO) plots a azimuth/range grid from the
%     lat/lon point given by EVLA/EVLO on the map with axes AX.  EVLA &
%     EVLO must be scalar and in units of degrees.  AX must be a handle for
%     a map produced with MMAP or MAPSTATIONS.  The grid steps at 10deg
%     increments in range and 15deg increments in azimuth.
%
%     MAPEVENTGRID(AX,EVLA,EVLO,C) sets the grid color.  C by default is
%     'b' (blue).
%
%     MAPEVENTGRID(AX,EVLA,EVLO,C,RRMAJ,RRMIN) defines the major & minor
%     range rings.  RRMAJ by default is [10 30:30:150 170].  RRMIN by
%     default is [20 40:30:130 50:30:140 160].
%
%     MAPEVENTGRID(AX,EVLA,EVLO,C,RRMAJ,RRMIN,AZMAJ,AZMIN) defines the
%     major & minor azimuth lines.  AZMAJ by default is [0:45:315].  AZMIN
%     by default is [15:45:330 30:45:345].
%
%    Notes:
%     - Azimuthal lines are tagged as 'az_minor' & 'az_major'.
%     - Range rings are tagged as 'rr_minor' & 'rr_major'.
%
%    Examples:
%     % Quickly make a map and place an azimuth/range grid on it:
%     mapeventgrid(mmap,0,0);
%
%     % This replaces the old MAPSTATIONS2:
%     ax=mapstations(data);
%     ev=getheader(data(1),'ev');
%     mapeventgrid(ax,ev(1),ev(2));
%
%    See also: MAPSTATIONS, MMAP, MAPFEATURE, MAPCMTS

%     Version History:
%        June 26, 2010 - initial version
%        Mar.  6, 2011 - changed line width for aesthetics
%        May  15, 2011 - updated docs 
%        Apr.  3, 2012 - minor doc update
%        Aug. 28, 2013 - do not require axes to be tagged
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 28, 2013 at 03:15 GMT

% todo:

% check nargin
error(nargchk(3,8,nargin));

% check ax
if(isempty(ax) || ~isscalar(ax) || ~isreal(ax) ...
        || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
    error('seizmo:mapeventgrid:badAxes',...
        'Map axes does not exist or is invalid!');
else
    axes(ax);
end

% check event position
if(~isreal(evla) || ~isscalar(evla))
    error('seizmo:mapeventgrid:badInput',...
        'EVLA must be a real-valued scalar!');
elseif(~isreal(evlo) || ~isscalar(evlo))
    error('seizmo:mapeventgrid:badInput',...
        'EVLO must be a real-valued scalar!');
end

% fix event position
[evla,evlo]=fixlatlon(evla,evlo);

% defaults
if(nargin<4 || isempty(c)); c='b'; end
if(nargin<5 || isempty(rrmajor)); rrmajor=[10 30:30:150 170]; end
if(nargin<6 || isempty(rrminor)); rrminor=[20 40:30:130 50:30:140 160]; end
if(nargin<7 || isempty(azmajor)); azmajor=0:45:315; end
if(nargin<8 || isempty(azminor)); azminor=[15:45:330 30:45:345]; end

% check optional args
if(~ischar(c) && (~isreal(c) || (isreal(c) && ~isequal(size(c),[1 3]))))
    error('seizmo:mapeventgrid:badInput',...
        'C must be a colorname or a RGB triplet!');
elseif(~isreal(rrmajor) || ~isreal(rrminor) ...
        || any(rrmajor<0 | rrmajor>180) || any(rrminor<0 | rrminor>180))
    error('seizmo:mapeventgrid:badInput',...
        'RRMAJOR & RRMINOR must be real-valued arrays within [0 180]!');
elseif(~isreal(azmajor) || ~isreal(azminor) ...
        || any(azmajor<0 | azmajor>360) || any(azminor<0 | azminor>360))
    error('seizmo:mapeventgrid:badInput',...
        'AZMAJOR & AZMINOR must be real-valued arrays within [0 360]!');
end

% number of points for lines
npts=200;

% get min/max range rings
rrmin=min([rrmajor(:); rrminor(:)]);
rrmax=max([rrmajor(:); rrminor(:)]);

% plot azimuthal lines
d2r=pi/180;
if(numel(azminor))
    [azsla,azslo]=sphericalfwd(evla,evlo,rrmin,azminor);
    [azela,azelo]=sphericalfwd(evla,evlo,rrmax,azminor);
    [azmla,azmlo]=sphericalfwd(evla,evlo,sum(rrmax+rrmin)/2,azminor);
    [lat1,lon1]=gcarc2latlon(azsla,azslo,azmla,azmlo,npts);
    [lat2,lon2]=gcarc2latlon(azela,azelo,azmla,azmlo,npts);
    lon=unwrap([lon1; lon2].*d2r,[],2)./d2r; % avoid streak from wraparound
    h=m_line(lon',[lat1; lat2]','color',c,'linewi',.25);
    set(h,'tag','az_minor');
end
if(numel(azmajor))
    [azsla,azslo]=sphericalfwd(evla,evlo,rrmin,azmajor);
    [azela,azelo]=sphericalfwd(evla,evlo,rrmax,azmajor);
    [azmla,azmlo]=sphericalfwd(evla,evlo,sum(rrmax+rrmin)/2,azmajor);
    [lat1,lon1]=gcarc2latlon(azsla,azslo,azmla,azmlo,npts);
    [lat2,lon2]=gcarc2latlon(azela,azelo,azmla,azmlo,npts);
    lon=unwrap([lon1; lon2].*d2r,[],2)./d2r; % avoid streak from wraparound
    h=m_line(lon',[lat1; lat2]','color',c,'linewi',.75);
    set(h,'tag','az_major');
end

% plot range rings
d2r=pi/180;
if(numel(rrminor))
    h=m_range_ring(evlo,evla,6371*d2r*rrminor,npts,'color',c,'linewi',.25);
    set(h,'tag','rr_minor');
end
if(numel(rrmajor))
    h=m_range_ring(evlo,evla,6371*d2r*rrmajor,npts,'color',c,'linewi',.75);
    set(h,'tag','rr_major');
end

end
