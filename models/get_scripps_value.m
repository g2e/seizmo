function [val]=get_scripps_value(s,lat,lon,depth)
%GET_SCRIPPS_VALUE    Returns dv value for a Scripps mantle model
%
%    Usage:    dv=get_scripps_value(mod,lat,lon,depth)
%
%    Description:
%     DV=GET_SCRIPPS_VALUE(MOD,LAT,LON,DEPTH) finds velocity values
%     associated with a Scripps mantle model MOD at the locations given by
%     LAT, LON, & DEPTH.  MOD should be created by SCRIPPS2MATLAB.  LAT,
%     LON, & DEPTH should be real-valued, equal valued arrays or scalars.
%     LAT & LON must be in degrees.  DEPTH must be in kilometers.
%
%    Notes:
%
%    Examples:
%     % Get a depth slice at 2800km (with an extra 90deg of longitude):
%     hmsl06p=load('HMSL06P');
%     [x,y]=meshgrid(0:2:450,90:-2:-90);
%     imagesc(0:2:450,90:-2:-90,get_scripps_value(hmsl06p,y,x,2800));
%     axis xy equal tight;
%     colormap(seis);
%
%    See also: SCRIPPS2MATLAB, READ_SCRIPPS_BINARY

%     Version History:
%        May  16, 2010 - initial version
%        June  1, 2010 - dv rather than dv%, fixed example
%        July  8, 2010 - really fixed example
%        Mar. 24, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 24, 2012 at 18:35 GMT

% todo

% check nargin
error(nargchk(4,4,nargin));

% check inputs
reqfields={'deplimits' 'latlimits' 'lonblks_at_lat' 'lat_block_idx' ...
    'lonlimits' 'model' 'blocks_per_layer' 'nlayers' 'blocksize'};
if(~isstruct(s) || any(~ismember(reqfields,fieldnames(s))))
    error('seizmo:get_scripps_value:badInput',...
        'S must be a Scripps struct created with SCRIPPS2MATLAB!');
elseif(~isreal(lat) || ~isreal(lon) || ~isreal(depth))
    error('seizmo:get_scripps_value:badInput',...
        'LAT, LON, & DEPTH must be real-valued arrays!');
elseif(~isequalsizeorscalar(lat,lon,depth))
    error('seizmo:get_scripps_value:badInput',...
        'LAT, LON, & DEPTH must be equal sized or scalar!');
end

% make all the same size
[lat,lon,depth]=expandscalars(lat,lon,depth);

% make sure lat/lon values are in proper ranges
[lat,lon]=fixlatlon(lat,lon);
lon(lon<0)=lon(lon<0)+360; % 0 to 360

% get sizes, make row vectors
sz=size(lat);
npts=numel(lat);
lat=lat(:).';
lon=lon(:).';
depth=depth(:).';

% get layer for each point
nd=s.nlayers;
[found,layer]=max(depth(ones(nd,1),:)<=s.deplimits(:,ones(npts,1)) ...
    & depth(ones(nd,1),:)>=s.deplimits(:,2*ones(npts,1)));

% get latitude strip for each point
nlat=180/s.blocksize;
[latidx,latidx]=max(lat(ones(nlat,1),:)<=s.latlimits(:,ones(npts,1)) ...
    & lat(ones(nlat,1),:)>=s.latlimits(:,2*ones(npts,1)));

% now get longitude index
nlon=s.lonblks_at_lat(latidx);
lonidx=nan(1,npts);
for i=1:npts
    [lonidx(i),lonidx(i)]=max(lon(ones(nlon(i),1),i)>=s.lonlimits(...
        s.lat_block_idx(latidx(i),1):s.lat_block_idx(latidx(i),2),1) ...
        & lon(ones(nlon(i),1),i)<=s.lonlimits(s.lat_block_idx(...
        latidx(i),1):s.lat_block_idx(latidx(i),2),2));
end

% put it all together
% (note found variable sets points outside range to 0)
val=reshape(found.'.*s.model((layer-1)*s.blocks_per_layer...
    +s.lat_block_idx(latidx,1).'+(lonidx-1)),sz);

end
