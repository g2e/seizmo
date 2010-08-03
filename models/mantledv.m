function [dlnv]=mantledv(model,lat,lon,depth)
%MANTLEDV    Returns the seismic velocity deviation for a mantle model
%
%    Usage:    dlnv=mantledv(model,lat,lon,depth)
%
%    Description: DLNV=MANTLEDV(MODEL,LAT,LON,DEPTH) returns the seismic
%     velocity deviation (from a 1D reference) for the mantle model MODEL
%     at the locations given by LAT/LON/DEPTH.  MODEL must be a string
%     matching one of the models output by AVAILABLE_3DMODELS.  LAT & LON
%     are in units of degrees.  DEPTH is in kilometers.  DLNV is the
%     fractional deviation from the reference velocity: dv/v0.
%
%    Notes:
%     - Model info is cached to speed up subsequent calls
%
%    Examples:
%     Compare several mantle models at 0deg lat, 0deg lon:
%      lat=0; lon=0;
%      depths=(50:50:2850)';
%      mod1=mantledv('sb4l18',lat,lon,depths);
%      mod2=mantledv('hmsl06s',lat,lon,depths);
%      mod3=mantledv('hmsl06p',lat,lon,depths);
%      figure;
%      plot([mod1 mod2 mod3],depths);
%      legend({'SB4L18' 'HMSL06S' 'HMSL06P'});
%      title(['dV @ LAT: ' num2str(lat) 'deg  LON: ' num2str(lon) 'deg']);
%
%    See also: AVAILABLE_3DMODELS, GET_SCRIPPS_VALUE

%     Version History:
%        June  1, 2010 - initial version
%        June  2, 2010 - changed variable name to reflect the truth,
%                        adjusts deplimits for variability in moho/cmb
%        July 25, 2010 - added several models
%        Aug.  2, 2010 - add mit-p08, minor dz04 fix, change cache local
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  2, 2010 at 16:05 GMT

% todo:
% - sph harm models
% - ani models

% check nargin
error(nargchk(4,4,nargin));

% use global for model caching
global SEIZMO

% check inputs
sm=size(model);
if(~ischar(model) || numel(sm)~=2 || sm(1)~=1)
    error('seizmo:mantledv:badInput',...
        'MODEL must be a string!');
elseif(~all(cellfun('isreal',{lat lon depth})))
    error('seizmo:mantledv:badInput',...
        'LAT/LON/DEPTH must be real-valued!');
elseif(~isequalsizeorscalar(lat,lon,depth))
    error('seizmo:mantledv:badInput',...
        'LAT/LON/DEPTH must be equal-sized or scalar!');
end

% expand scalar inputs
[lat,lon,depth]=expandscalars(lat,lon,depth);

% make sure lat/lon values are in proper ranges
[lat,lon]=fixlatlon(lat,lon);

% act based on mantle model
switch lower(model)
    case {'dz04'}
        % load cached model if it exists
        try
            model=SEIZMO.MANTLEMODEL.DZ04;
        catch
            % not there so load and cache
            model=load('DZ04');
            model.dvp=model.dvp/100;
            SEIZMO.MANTLEMODEL.DZ04=model;
        end
        
        % get dlnv
        dlnv=get_interp_value(model,lat,lon,depth);
    case {'sb4l18'}
        % load cached model if it exists
        try
            model=SEIZMO.MANTLEMODEL.SB4L18;
        catch
            % not there so load and cache
            model=load('SB4L18');
            SEIZMO.MANTLEMODEL.SB4L18=model;
        end
        
        % allow depths between 0 & 2900 to allow for variable moho and cmb
        model.deplimits(end)=0;
        model.deplimits(1)=2900;
        
        % get dlnv
        dlnv=get_scripps_value(model,lat,lon,depth);
    case {'hmsl06p' 'hmslp06' 'hmsl-p06' 'hmsl-06p'}
        % load cached model if it exists
        try
            model=SEIZMO.MANTLEMODEL.HMSL06P;
        catch
            % not there so load and cache
            model=load('HMSL06P');
            SEIZMO.MANTLEMODEL.HMSL06P=model;
        end
        
        % allow depths between 0 & 2900 to allow for variable moho and cmb
        model.deplimits(end)=0;
        model.deplimits(1)=2900;
        
        % get dlnv
        dlnv=get_scripps_value(model,lat,lon,depth);
    case {'hmsl06s' 'hmsls06' 'hmsl-s06' 'hmsl-06s'}
        % load cached model if it exists
        try
            model=SEIZMO.MANTLEMODEL.HMSL06S;
        catch
            % not there so load and cache
            model=load('HMSL06S');
            SEIZMO.MANTLEMODEL.HMSL06S=model;
        end
        
        % allow depths between 0 & 2900 to allow for variable moho and cmb
        model.deplimits(end)=0;
        model.deplimits(1)=2900;
        
        % get dlnv
        dlnv=get_scripps_value(model,lat,lon,depth);
    case {'mit08' 'mit08p' 'mit-p08' 'mitp08'}
        % load cached model if it exists
        try
            model=SEIZMO.MANTLEMODEL.MIT08;
        catch
            % not there so load and cache
            model=load('MIT08');
            model.dvp=model.dvp/100;
            SEIZMO.MANTLEMODEL.MIT08=model;
        end
        
        % get dlnv
        dlnv=get_mit_value(model,lat,lon,depth);
    case {'pri05p' 'pri-p05' 'pri-05p' 'prip05'}
        % load cached model if it exists
        try
            model=SEIZMO.MANTLEMODEL.PRI05P;
        catch
            % not there so load and cache
            model=load('PRI05','ppts','refmod','name','reference');
            model.name=model.name.dvp;
            SEIZMO.MANTLEMODEL.PRI05P=model;
        end
        
        % get dlnv
        dlnv=get_princeton_value(model,lat,lon,depth);
    case {'pri05s' 'pri-s05' 'pri-05s' 'pris05'}
        % load cached model if it exists
        try
            model=SEIZMO.MANTLEMODEL.PRI05S;
        catch
            % not there so load and cache
            model=load('PRI05','spts','refmod','name','reference');
            model.name=model.name.dvs;
            SEIZMO.MANTLEMODEL.PRI05S=model;
        end
        
        % get dlnv
        dlnv=get_princeton_value(model,lat,lon,depth);
    case {'s20rtsb' 's20rts'}
        % load cached model if it exists
        try
            model=SEIZMO.MANTLEMODEL.S20RTS;
        catch
            % not there so load and cache
            model=load('S20RTS');
            model.dvs=model.dvs/100;
            SEIZMO.MANTLEMODEL.S20RTS=model;
        end
        
        % get dlnv
        dlnv=get_interp_value(model,lat,lon,depth);
    case {'saw24b16'}
        % load cached model if it exists
        try
            model=SEIZMO.MANTLEMODEL.SAW24B16;
        catch
            % not there so load and cache
            model=load('SAW24B16');
            model.dvs=model.dvs/100;
            SEIZMO.MANTLEMODEL.SAW24B16=model;
        end
        
        % get dlnv
        dlnv=get_interp_value(model,lat,lon,depth);
    case {'tx2006'}
        % load cached model if it exists
        try
            model=SEIZMO.MANTLEMODEL.TX2006;
        catch
            % not there so load and cache
            model=load('TX2006');
            model.dvs=model.dvs/100;
            SEIZMO.MANTLEMODEL.TX2006=model;
        end
        
        % get dlnv
        dlnv=get_texas_value(model,lat,lon,depth);
    case {'tx2007'}
        % load cached model if it exists
        try
            model=SEIZMO.MANTLEMODEL.TX2007;
        catch
            % not there so load and cache
            model=load('TX2007');
            model.dvs=model.dvs/100;
            SEIZMO.MANTLEMODEL.TX2007=model;
        end
        
        % get dlnv
        dlnv=get_texas_value(model,lat,lon,depth);
    otherwise
        error('seizmo:mantledv:badModelName',...
            'Unknown mantle model: %s',model);
end

end

function dlnv=get_texas_value(model,lat,lon,depth)
% get dlnv
% - find layer first (b/c they are of variable thickness)
% - nearest neighbor interpolation
dlnv=nan(size(depth));
for i=1:size(model.depth,1)
    idx=depth>=model.depth(i,1) & depth<=model.depth(i,2);
    if(any(idx(:)))
        dlnv(idx)=interp2(model.lon,model.lat,model.dvs(:,:,i),...
            lon(idx),lat(idx),'*nearest');
    end
end
end

function dlnv=get_interp_value(model,lat,lon,depth)
% get dlnv
% - spline interpolation
if(isfield(model,'dvs'))
    dlnv=interp3(model.lon,model.lat,model.depth,model.dvs,...
        lon,lat,depth,'spline');
else % dvp
    dlnv=interp3(model.lon,model.lat,model.depth,model.dvp,...
        lon,lat,depth,'spline');
end
end

function dlnv=get_princeton_value(model,lat,lon,depth)
% get dlnv
% - natural interpolation from tetrahedral mesh
if(isfield(model,'spts'))
    dt=DelaunayTri(model.spts(:,1:3));
    f=TriScatteredInterp(dt,model.spts(:,4),'natural');
else % ppts
    dt=DelaunayTri(model.ppts(:,1:3));
    f=TriScatteredInterp(dt,model.ppts(:,4),'natural');
end
[x,y,z]=geographic2xyz(lat,lon,depth);
dlnv=f(x,y,z);
end

function dlnv=get_mit_value(model,lat,lon,depth)

% fix lon
lon(lon<0)=lon(lon<0)+360;

% get dlnv
% - nearest neighbor interpolation
% - get linear indices
lat=floor((90-lat)/(180/256))+1;
lat(lat==257)=256; % handle south pole
lon=floor(lon/(360/512));
lon(lon==512)=511; % handle meridian
depth=floor(depth/(2891.5/64));
depth(depth<0)=0; % handle surface
depth(depth>63)=63; % handle cmb
dlnv=model.dvp(lat+lon*256+depth*256*512);

% this would be ok if the south pole returned non-nan ...
%dlnv=interp3(model.lon,model.lat,model.depth,model.dvp,...
%    lon,lat,depth,'*nearest');

end
