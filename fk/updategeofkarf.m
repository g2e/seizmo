function [varargout]=updategeofkarf(arf,ax)
%UPDATEGEOFKARF    Quickly updates an existing geofkarf plot with a new arf
%
%    Usage:    updategeofkarf(arf,ax)
%              ax=updategeofkarf(arf)
%
%    Description:
%     UPDATEGEOFKARF(ARF,AX) draws a new geofkarf map given by ARF in an
%     existing axes AX.  This is mainly intended for exploring geofkarf
%     datasets and for making movies in a faster fashion.
%
%     AX=UPDATEGEOFKARF(ARF) is the same as calling PLOTGEOFKARF(ARF) --
%     ie. a new figure is drawn.
%
%    Notes:
%
%    Examples:
%
%    See also: GEOFKARF, PLOTGEOFKARF, GEOFKARF2MAP, GEOFKSUBARF,
%              GEOFKARFSLOWSLIDE, CHKGEOFKARFSTRUCT

%     Version History:
%        July  8, 2010 - update for new struct
%        Oct.  6, 2010 - truncate title if too many ARF locations
%        Dec.  8, 2010 - use '^o' for deg symbol rather than \circ
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 10:45 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check fk struct
error(chkgeofkarfstruct(arf));

% don't allow array/volume
if(~isscalar(arf) || any(arf.volume))
    error('seizmo:updategeofkarf:badInput',...
        'ARF must be a scalar geofkarf struct and not a volume!');
end

% just replot if ax isn't an axes handle
if(nargin<2 || ~isscalar(ax) || ~isreal(ax) || ~ishandle(ax) ...
        || ~strcmp('axes',get(ax,'type')))
    ax=plotgeofkarf(arf);
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

% rescale arf
switch zerodb
    case 'min'
        arf.beam=arf.beam-min(arf.beam(:));
        arf.normdb=arf.normdb+min(arf.beam(:));
    case 'max'
        arf.beam=arf.beam-max(arf.beam(:));
        arf.normdb=arf.normdb+max(arf.beam(:));
    case 'median'
        arf.beam=arf.beam-median(arf.beam(:));
        arf.normdb=arf.normdb+median(arf.beam(:));
    case 'abs'
        arf.beam=arf.beam+arf.normdb;
        arf.normdb=0;
end

% reshape arf & account for pcolor
nlat=numel(unique(arf.latlon(:,1)));
nlon=numel(unique(arf.latlon(:,2)));
arf.latlon=reshape(arf.latlon,[nlon nlat 2]);
arf.beam=reshape(arf.beam,[nlon nlat]);
latstep=arf.latlon(1,2,1)-arf.latlon(1,1,1);
lonstep=arf.latlon(2,1,2)-arf.latlon(1,1,2);
arf.latlon(:,:,1)=arf.latlon(:,:,1)-latstep/2;
arf.latlon(:,:,2)=arf.latlon(:,:,2)-lonstep/2;
arf.latlon=arf.latlon([1:end end],[1:end end],:);
arf.latlon(:,end,1)=arf.latlon(:,end,1)+latstep;
arf.latlon(end,:,2)=arf.latlon(end,:,2)+lonstep;
arf.beam=arf.beam([1:end end],[1:end end]);

% convert to map coordinates
[arf.latlon(:,:,2),arf.latlon(:,:,1)]=m_ll2xy(...
    arf.latlon(:,:,2),arf.latlon(:,:,1),'clip','patch');

% find previous
pc=findobj(ax,'tag','m_pcolor');

% slip in new data
set(pc(1),...
    'xdata',arf.latlon(:,:,2),'ydata',arf.latlon(:,:,1),...
    'zdata',0*arf.latlon(:,:,2),'cdata',double(arf.beam));

% adjust title
smn=min(arf.horzslow); smx=max(arf.horzslow);
if(arf.nsw<=5)
    titstr=cell(arf.nsw,1);
    for i=1:arf.nsw
        titstr{i}=sprintf(['SLOWNESS: %gs/^o, LAT: %g^o, ' ...
            'LON: %g^o, PERIOD: %gs'],arf.horzslow0(i),...
            arf.latlon0(i,1),arf.latlon0(i,2),1/arf.freq0(i));
    end
else
    titstr{1}=[num2str(arf.nsw) ' Locations'];
end
set(get(ax,'Title'),'string',...
    [{[]}; 'Array Response Function @ '; titstr; ...
    ['Number of Stations: ' num2str(arf.nsta)]; ...
    ['Horiz. Slowness : ' num2str(smn) ' to ' num2str(smx) ' s/^o']; ...
    ['0 dB = ' num2str(arf.normdb) 'dB']; {[]}]);

if(nargout); varargout{1}=ax; end

end
