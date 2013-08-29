function []=updategeofkmap(map,ax)
%UPDATEGEOFKMAP    Quickly updates an existing geofk plot with a new map
%
%    Usage:    updategeofkmap(map,ax)
%
%    Description:
%     UPDATEGEOFKMAP(MAP,AX) quickly updates an existing geofk map given
%     with axes handle AX with the new beam data in MAP.  This is intended
%     for exploring geofk datasets by making "movies" and is used by
%     functions such as GEOFKSLOWSLIDE and GEOFKFREQSLIDE.
%
%    Notes:
%
%    Examples:
%     % Also useful if you just want to switch between two datasets:
%     map(1)=geofkvol2map(geofkxcvolume(data,[la(:) lo(:)],30,[0.1 0.2]));
%     map(2)=geofkvol2map(geofkxcvolume(data,[la(:) lo(:)],27,[0.1 0.2]));
%     ax=plotgeofkmap(map(1));
%     for i=2:100
%         pause(.25);
%         updategeofkmap(map(2-mod(i,2)),ax);
%     end
%
%    See also: PLOTGEOFKMAP, GEOFKXCVOLUME, GEOFKFREQSLIDE, GEOFKSLOWSLIDE,
%              GEOFKXCHORZVOLUME, CHKGEOFKSTRUCT

%     Version History:
%        June 25, 2010 - initial version
%        July  6, 2010 - update for new struct
%        Dec.  8, 2010 - improved degrees symbol usage (^o)
%        Apr.  4, 2012 - minor doc update
%        Aug. 28, 2013 - require axes, error if axes not valid, better
%                        checking, use mmap to update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 28, 2013 at 19:05 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check fk struct
error(chkgeofkstruct(map));

% don't allow array/volume
if(~isscalar(map) || any(map.volume))
    error('seizmo:updategeofkmap:badInput',...
        'MAP must be a scalar geofk struct and not a volume!');
end

% just replot if ax isn't an axes handle
if(nargin<2 || ~isscalar(ax) || ~isreal(ax) || ~ishandle(ax) ...
        || ~strcmp('axes',get(ax,'type')) ...
        || ~strcmp('geofkmap',get(ax,'tag')))
    error('seizmo:updategeofkmap:badAxes',...
        'Map axes for updating geofk spectra does not exist!');
end

% get zerodb/dblim
userdata=get(ax,'userdata');

% rescale beam
switch userdata.zerodb
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

% find previous objects & delete
delete(findobj(ax,'tag','m_pcolor'));
delete(findobj(ax,'tag','stations'));

% update map
held=ishold(ax);
if(~held); hold(ax,'on'); end
mmap('image',{map.latlon(:,:,1) map.latlon(:,:,2) double(map.beam)},...
    'st',[map.stla map.stlo],'parent',ax);
if(~held); hold(ax,'off'); end

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

end
