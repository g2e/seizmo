function [crust]=getcrust(lat,lon,model)
%GETCRUST    Returns crustal info at specified location(s)
%
%    Usage:    crust=getcrust(lat,lon)
%              crust=getcrust(lat,lon,model)
%
%    Description:
%     CRUST=GETCRUST(LAT,LON) returns struct CRUST containing info about
%     the crustal structure at the positions given by LAT & LON.  LAT & LON
%     are expected to be scalars or equal sized arrays and have units in
%     degrees.  CRUST has a layout as follows:
%      CRUST.model       -- crustal model name
%           .top         -- elev from sealevel to layer top in km (Nx9)
%           .vp          -- isotropic p-velocities in km/s (Nx9)
%           .vs          -- isotropic s-velocities in km/s (Nx9)
%           .rho         -- densities in kg/m3 (Nx9)
%           .thk         -- layer kilometer thicknesses (Nx8)
%     where N is the number of LAT/LON points.  The values come from the
%     Crust1.0 model which is decribed in more detail in the Notes section.
%
%     CRUST=GETCRUST(LAT,LON,MODEL) allows specifying a different crustal
%     model.  Currently available models:
%      'CRUST1.0', 'CRUST2.0', 'CRUST5.1'
%     Note that CRUST2.0 & CRUST5.1 have one less sedimentary layers
%     compared to CRUST1.0.  This missing layer is assigned as a 0km
%     thickness middle sedimentary layer (see Notes section).
%
%    Notes:
%     - Layers are sorted as follows by columns:
%
%                             SEDS                 CRUST
%       WATER     ICE  UPPER MIDDLE  LOWER  UPPER MIDDLE  LOWER   MOHO
%       top(1) top(2) top(3) top(4) top(5) top(6) top(7) top(8) top(9)
%       vp(1)   vp(2)  vp(3)  vp(4)  vp(5)  vp(6)  vp(7)  vp(8)  vp(9)
%       vs(1)   vs(2)  vs(3)  vs(4)  vs(5)  vs(6)  vs(7)  vs(8)  vs(9)
%       rho(1) rho(2) rho(3) rho(4) rho(5) rho(6) rho(7) rho(8) rho(9)
%       thk(1) thk(2) thk(3) thk(4) thk(5) thk(6) thk(7) thk(8) ------
%
%    Examples:
%     % What is the crustal structure below Saint Louis?
%     getcrust(38.649,-90.305)
%
%     % Moho map:
%     [lon,lat]=meshgrid(-179.5:179.5,-89.5:89.5);
%     c=getcrust(lat,lon);
%     fh=figure;
%     ax=axes('parent',fh);
%     imagesc(-179.5:179.5,-89.5:89.5,...
%         reshape(-c.top(:,9),[180 360]),'parent',ax);
%     axis(ax,'xy');
%     cbh=colorbar('peer',ax);
%     ylabel(cbh,'km');
%     title(ax,'Moho Depth');
%
%     % Elevation & Bathymetry map:
%     [lon,lat]=meshgrid(-179.5:179.5,-89.5:89.5);
%     c=getcrust(lat,lon);
%     fh=figure;
%     ax=axes('parent',fh);
%     imagesc(-179.5:179.5,-89.5:89.5,...
%         reshape(c.top(:,2)*1000,[180 360]),'parent',ax);
%     axis(ax,'xy');
%     cbh=colorbar('peer',ax);
%     colormap(ax,topo_colormap(c.top(:,2)*1000));
%     ylabel(cbh,'meters');
%     title(ax,'Elevation & Bathymetry');
%
%     % Crustal Thickness map:
%     [lon,lat]=meshgrid(-179.5:179.5,-89.5:89.5);
%     c=getcrust(lat,lon);
%     fh=figure;
%     ax=axes('parent',fh);
%     imagesc(-179.5:179.5,-89.5:89.5,...
%         reshape(sum(c.thk(:,2:8),2),[180 360]),'parent',ax);
%     axis(ax,'xy');
%     cbh=colorbar('peer',ax);
%     ylabel(cbh,'km');
%     title(ax,'Crustal Thickness (w/ Ice)');
%
%    See also: CRUCOR, MANCOR, CRUSTLESS_RAYPATHS

%     Version History:
%        Jan. 23, 2014 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 23, 2014 at 18:35 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% use global for model caching
global SEIZMO

% check lat/lon
if(~isreal(lat) || ~isreal(lon))
    error('seizmo:getcrust:badInput',...
        'LAT & LON must be real-valued arrays!');
elseif(~isequalsizeorscalar(lat,lon))
    error('seizmo:getcrust:badInput',...
        'LAT & LON must be equal sized or scalar!');
end

% make the same size
[lat,lon]=expandscalars(lat,lon);

% make sure lat/lon values are in proper ranges
[lat,lon]=fixlatlon(lat,lon);

% make column vectors
lat=lat(:);
lon=lon(:);

% default/check model
if(nargin<3 || isempty(model)); model='CRUST1.0'; end
if(~ischar(model) || size(model,1)~=1 || ndims(model)~=2)
    error('seizmo:getcrust:badInput',...
        'MODEL must be a string!');
end
switch lower(model)
    case {'crust1.0' 'c1' 'c1.0' '1.0' '1'}
        % load cached model if it exists
        try
            mod=SEIZMO.CRUSTALMODEL.CRUST10I;
        catch
            % not there so load and cache
            mod=load('CRUST10I');
            SEIZMO.CRUSTALMODEL.CRUST10I=mod;
        end
        crust.model='CRUST1.0';
        nlat=180;
        nlon=360;
        step=1;
    case {'crust2.0' 'c2' 'c2.0' '2.0' '2'}
        % load cached model if it exists
        try
            mod=SEIZMO.CRUSTALMODEL.CRUST20I;
        catch
            % not there so load and cache
            mod=load('CRUST20I');
            SEIZMO.CRUSTALMODEL.CRUST20I=mod;
        end
        crust.model='CRUST2.0';
        nlat=90;
        nlon=180;
        step=2;
    case {'crust5.1' 'c5' 'c5.1' '5.1' '5'}
        % load cached model if it exists
        try
            mod=SEIZMO.CRUSTALMODEL.CRUST51I;
        catch
            % not there so load and cache
            mod=load('CRUST51I');
            SEIZMO.CRUSTALMODEL.CRUST51I=mod;
        end
        crust.model='CRUST5.1';
        nlat=36;
        nlon=72;
        step=5;
    otherwise
        error('seizmo:getcrust:badInput',...
            'MODEL is unknown!');
end

% lat/lon to i/j
i=fix((90-lat)/step)+1;
i(i==(nlat+1))=nlat;
j=fix((180+lon)/step)+1;
j(j==(nlon+1))=nlon;
idx=(j-1)*nlat+i;

% get crust
crust.top=double(mod.top(idx,:))/100;
crust.vp=double(mod.vp(idx,:))/100;
crust.vs=double(mod.vs(idx,:))/100;
crust.rho=double(mod.rho(idx,:))/100;
crust.thk=-diff(crust.top,1,2);

end
