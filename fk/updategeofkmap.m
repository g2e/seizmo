function [varargout]=plotgeofssupdate(s,ax)
%PLOTGEOFSSUPDATE    Updates an existing geofss plot with a new spectra
%
%    Usage:    plotgeofssupdate(s,ax)
%              ax=plotgeofssupdate(s)
%
%    Description:
%     UPDATEGEOFSS(S,AX) draws a new geofk beam data map given by MAP
%     in an existing axes AX.  This is mainly intended for exploring geofk
%     datasets and for making movies in a faster fashion.
%
%     AX=UPDATEGEOFKMAP(MAP) is the same as calling PLOTGEOFKMAP(MAP) --
%     ie. a new figure is drawn.
%
%    Notes:
%
%    Examples:
%
%    See also: PLOTGEOFKMAP, GEOFKXCVOLUME, GEOFKFREQSLIDE, GEOFKSLOWSLIDE,
%              GEOFKXCHORZVOLUME, CHKGEOFKSTRUCT

%     Version History:
%        June 25, 2010 - initial version
%        July  6, 2010 - update for new struct
%        Dec.  8, 2010 - improved degrees symbol usage (^o)
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 19:05 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check fk struct
error(chkgeofkstruct(map));

% don't allow array/volume
if(~isscalar(map) || any(map.volume))
    error('seizmo:updategeofkmap:badInput',...
        'MAP must be a scalar geofk struct and not a volume!');
end

% just replot if ax isn't an axes handle
if(nargin<2 || ~isscalar(ax) || ~isreal(ax) || ~ishandle(ax) ...
        || ~strcmp('axes',get(ax,'type')))
    ax=plotgeofkmap(map);
    if(nargout); varargout{1}=ax; end
    return;
end

% get zerodb/dblim
userdata=get(ax,'userdata');
if(isempty(userdata) || ~isstruct(userdata) ...
        || any(~isfield(userdata,{'zerodb' 'dblim'})))
    %dblim=[-12 0];
    zerodb='max';
else
    %dblim=userdata.dblim;
    zerodb=userdata.zerodb;
end

% rescale beam
switch zerodb
    case 'min'
        map.beam=map.beam-min(map.beam(:));
        map.normdb=map.normdb+min(map.beam(:));
    case 'max'
        map.beam=map.beam-max(map.beam(:));
        map.normdb=map.normdb+max(map.beam(:));
    case 'median'
        map.beam=map.beam-median(map.beam(:));
        map.normdb=map.normdb+median(map.beam(:));
    case 'abs'
        map.beam=map.beam+map.normdb;
        map.normdb=0;
end

% reshape beam data & account for pcolor
nlat=numel(unique(map.latlon(:,1)));
nlon=numel(unique(map.latlon(:,2)));
map.latlon=reshape(map.latlon,[nlon nlat 2]);
map.beam=reshape(map.beam,[nlon nlat]);
latstep=map.latlon(1,2,1)-map.latlon(1,1,1);
lonstep=map.latlon(2,1,2)-map.latlon(1,1,2);
map.latlon(:,:,1)=map.latlon(:,:,1)-latstep/2;
map.latlon(:,:,2)=map.latlon(:,:,2)-lonstep/2;
map.latlon=map.latlon([1:end end],[1:end end],:);
map.latlon(:,end,1)=map.latlon(:,end,1)+latstep;
map.latlon(end,:,2)=map.latlon(end,:,2)+lonstep;
map.beam=map.beam([1:end end],[1:end end]);

% convert to map coordinates
[map.latlon(:,:,2),map.latlon(:,:,1)]=m_ll2xy(...
    map.latlon(:,:,2),map.latlon(:,:,1),'clip','patch');

% find previous
pc=findobj(ax,'tag','m_pcolor');

% slip in new data
set(pc(1),...
    'xdata',map.latlon(:,:,2),'ydata',map.latlon(:,:,1),...
    'zdata',0*map.latlon(:,:,2),'cdata',double(map.beam));

% adjust title
fmin=min(map.freq); fmax=max(map.freq);
smn=min(map.horzslow); smx=max(map.horzslow);
set(get(ax,'Title'),'string',...
    {[] ['Number of Stations:  ' num2str(map.nsta)] ...
    ['Begin Time:  ' sprintf('%d.%03d %02d:%02d:%02g',map.butc) ' UTC'] ...
    ['End Time  :  ' sprintf('%d.%03d %02d:%02d:%02g',map.eutc) ' UTC'] ...
    ['Period    :  ' num2str(1/fmax) ' to ' num2str(1/fmin) ' s'] ...
    ['Horiz. Slowness :  ' num2str(smn) ' to ' num2str(smx) ' s/^o'] ...
    ['0 dB = ' num2str(map.normdb) 'dB'] []});

if(nargout); varargout{1}=ax; end

end
