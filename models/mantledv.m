function [dlnv]=mantledv(model,lat,lon,depth)
%MANTLEDV    Returns the seismic velocity deviation for a mantle model
%
%    Usage:    dlnv=mantledv(model,lat,lon,depth)
%
%    Description: DLNV=MANTLEDV(MODEL,LAT,LON,DEPTH) returns the seismic
%     velocity deviation (from a 1D reference) for the mantle model MODEL
%     at the locations given by LAT/LON/DEPTH.  MODEL must be a string
%     matching one of the models output by AVAILABLE_3DMODELS.  LAT & LON
%     are in units of degrees.  DEPTH is in kilometers.  DLNV is is the
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  2, 2010 at 16:05 GMT

% todo:

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

% act based on mantle model
switch lower(model)
    case {'sb4l18'}
        % load cached model if it exists
        try
            model=SEIZMO.MANTLEDV.SB4L18;
        catch
            % not there so load and cache
            model=load('SB4L18');
            SEIZMO.MANTLEDV.SB4L18=model;
        end
        
        % allow depths between 0 & 2900 to allow for variable moho and cmb
        model.deplimits(end)=0;
        model.deplimits(1)=2900;
        
        % get dlnv
        dlnv=get_scripps_value(model,lat,lon,depth);
    case {'hmsl06p' 'hmslp06' 'hmsl-p06' 'hmsl-06p'}
        % load cached model if it exists
        try
            model=SEIZMO.MANTLEDV.HMSL06P;
        catch
            % not there so load and cache
            model=load('HMSL06P');
            SEIZMO.MANTLEDV.HMSL06P=model;
        end
        
        % allow depths between 0 & 2900 to allow for variable moho and cmb
        model.deplimits(end)=0;
        model.deplimits(1)=2900;
        
        % get dlnv
        dlnv=get_scripps_value(model,lat,lon,depth);
    case {'hmsl06s' 'hmsls06' 'hmsl-s06' 'hmsl-06s'}
        % load cached model if it exists
        try
            model=SEIZMO.MANTLEDV.HMSL06S;
        catch
            % not there so load and cache
            model=load('HMSL06S');
            SEIZMO.MANTLEDV.HMSL06S=model;
        end
        
        % allow depths between 0 & 2900 to allow for variable moho and cmb
        model.deplimits(end)=0;
        model.deplimits(1)=2900;
        
        % get dlnv
        dlnv=get_scripps_value(model,lat,lon,depth);
    otherwise
        error('seizmo:mantledv:badModelName',...
            'Unknown mantle model: %s',model);
end

end
