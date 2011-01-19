function [varargout]=map_cmb_profiles(pf,field,clim,varargin)
%MAP_CMB_PROFILES    Map CMB profile measurements
%
%    Usage:    map_cmb_profiles(pf,field)
%              map_cmb_profiles(pf,field,clim)
%              map_cmb_profiles(pf,field,clim,...,'option',value,...)
%              ax=map_cmb_profiles(...)
%
%    Description:
%     MAP_CMB_PROFILES(PF,FIELD) plots the profiles in struct PF from
%     SLOWDECAYPAIRS or SLOWDECAYPROFILES, colored by the values given in
%     the field FIELD.  FIELD may be 'SLOW', 'CSLOW', 'DECAY' or 'CDECAY'.
%
%     MAP_CMB_PROFILES(PF,FIELD,CLIM) adjusts the color limits for the
%     profile coloring.  The default is [] (empty matrix) and will use the
%     data limits.  CLIM should be specified as [MIN MAX].
%
%     MAP_CMB_PROFILES(PF,FIELD,CLIM,...,'OPTION',VALUE,...) allows passing
%     options for mapping.  See MAPLOCATIONS for details.
%
%     AX=MAP_CMB_PROFILES(...) returns the axis handle AX of the plot.
%
%    Notes:
%
%    Examples:
%     % Plot corrected profile slowness with color limits of 4.4-4.7:
%     map_cmb_profiles(pf,'cslow',[4.4 4.7])
%
%    See also: MAPLOCATIONS, SLOWDECAYPAIRS, SLOWDECAYPROFILES,
%              CMB_1ST_PASS, CMB_OUTLIERS, CMB_2ND_PASS

%     Version History:
%        Dec. 12, 2010 - initial version
%        Dec. 13, 2010 - fix wrap-around issue
%        Jan. 18, 2011 - new pf fields, fix single profile bug
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 18, 2011 at 13:35 GMT

% todo:
% - handle single profiles

% check nargin
error(nargchk(2,inf,nargin));
if(nargin>3 && ~mod(nargin,2))
    error('seizmo:map_cmb_profiles:badNumInputs',...
        'Unpaired Option/Value!');
end

% check profiles struct
reqfields={'gcdist','azwidth','slow','slowerr','decay','decayerr',...
    'cslow','cslowerr','cdecay','cdecayerr','cluster','kname','st','ev',...
    'delaz','corrections','corrcoef','freq','phase','runname','dirname',...
    'time'};
if(~isstruct(pf) || any(~isfield(pf,reqfields)))
    error('seizmo:map_cmb_profiles:badInput',...
        ['PF must be a struct with the fields:\n' ...
        sprintf('''%s'' ',reqfields{:}) '!']);
end

% check field
if(~ischar(field) || ...
        ~any(strcmpi(field,{'slow' 'cslow' 'decay' 'cdecay'})))
    error('seizmo:map_cmb_profiles:badInput',...
        'FIELD must be ''SLOW'', ''CSLOW'', ''DECAY'', or ''CDECAY''!');
end
field=lower(field);

% slowness to color
slow=[pf.(field)].';
switch field
    case {'decay' 'cdecay'}
        % flip polarity for ease
        slow=-slow;
end
if(nargin<3 || isempty(clim))
    if(numel(pf)>1)
        clim=[min(slow) max(slow)];
    else
        % only fails if zero
        clim=[slow*.9 slow*1.1];
    end
end
sv=(slow-clim(1))/(clim(2)-clim(1));
sv(sv<0)=0;
sv(sv>1)=1;
color=[2*sv 2*sv 2-2*sv];
color(color(:,2)>1,2)=2-color(color(:,2)>1,2);
color(color>1)=1;

% unique stations & events
st=unique(cell2mat({pf.st}'),'rows');
ev=unique(cell2mat({pf.ev}'),'rows');

% call maplocations
ax=maplocations('stations',st(:,1:2),'events',ev(:,1:2),varargin{:});

% loop over profiles, plotting each
hold(ax,'on');
for a=1:numel(pf)
    % get profile endpoints
    evla=pf(a).ev(1,1);
    evlo=pf(a).ev(1,2);
    distmin=min(pf(a).delaz(:,4))-44*111.19;
    distmax=max(pf(a).delaz(:,4))-44*111.19;
    azavg=mean(pf(a).delaz(:,2));
    [lat1,lon1]=vincentyfwd(evla,evlo,distmin,azavg);
    [lat2,lon2]=vincentyfwd(evla,evlo,distmax,azavg);
    
    % plot profile
    [lat,lon]=gcarc2latlon(lat1,lon1,lat2,lon2);
    lon=unwrap(lon*pi/180,[],2)*180/pi; % avoid wraparound streak
    m_line(lon',lat','linewi',2,'color',color(a,:));
end
hold(ax,'off');

% add colormap
colormap(ax,blue2red);
set(ax,'clim',clim);
colorbar('peer',ax);

% move stations/events forward
sta=findobj(ax,'tag','stations');
evt=findobj(ax,'tag','events');
movekids(evt,'front');
movekids(sta,'front');

% optional output
if(nargout); varargout{1}=ax; end

end
