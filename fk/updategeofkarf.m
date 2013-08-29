function []=updategeofkarf(arf,ax)
%UPDATEGEOFKARF    Quickly updates an existing geofkarf plot with a new arf
%
%    Usage:    updategeofkarf(arf,ax)
%
%    Description:
%     UPDATEGEOFKARF(ARF,AX) quickly updates an existing geofkarf map given
%     with axes handle AX with the new beam response data in ARF.  This is
%     intended for exploring geofkarf datasets by making "movies" and is
%     used by functions such as GEOFKARFSLOWSLIDE.
%
%    Notes:
%
%    Examples:
%     % See what happens to the ARF as you delete 1 station from a global
%     % array (looping over all possible stations):
%     st=randlatlon(30);
%     ev=randlatlon(6);
%     [la,lo]=meshgrid(-90:90,-180:180);
%     arf=geofkarf(st,[la(:) lo(:)],30,ev,30,1/500,'center');
%     ax=plotgeofkarf(geofkarf2map(arf),[-6 0]);
%     for i=1:30
%         drawnow;
%         st1=st; st1(i,:)=[];
%         arf=geofkarf(st1,[la(:) lo(:)],30,ev,30,1/500,'center');
%         updategeofkarf(geofkarf2map(arf),ax);
%     end
%
%    See also: GEOFKARF, PLOTGEOFKARF, GEOFKARF2MAP, GEOFKSUBARF,
%              GEOFKARFSLOWSLIDE, CHKGEOFKARFSTRUCT

%     Version History:
%        July  8, 2010 - update for new struct
%        Oct.  6, 2010 - truncate title if too many ARF locations
%        Dec.  8, 2010 - use '^o' for deg symbol rather than \circ
%        Apr.  4, 2012 - minor doc update
%        Aug. 28, 2013 - require axes, error if axes not valid, better
%                        checking, use mmap to update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 28, 2013 at 10:45 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% check fk struct
error(chkgeofkarfstruct(arf));

% don't allow array/volume
if(~isscalar(arf) || any(arf.volume))
    error('seizmo:updategeofkarf:badInput',...
        'ARF must be a scalar geofkarf struct and not a volume!');
end

% just replot if ax isn't an axes handle
if(nargin<2 || ~isscalar(ax) || ~isreal(ax) || ~ishandle(ax) ...
        || ~strcmp('axes',get(ax,'type')) ...
        || ~strcmp('geofkarf',get(ax,'tag')))
    error('seizmo:updategeofkarf:badAxes',...
        'Map axes for updating geofkarf spectra does not exist!');
end

% get zerodb/dblim
userdata=get(ax,'userdata');

% rescale arf
switch userdata.zerodb
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

% reshape beam
nlat=numel(unique(arf.latlon(:,1)));
nlon=numel(unique(arf.latlon(:,2)));
arf.latlon=reshape(arf.latlon,[nlon nlat 2]);
arf.beam=reshape(arf.beam,[nlon nlat]);

% find previous objects & delete
delete(findobj(ax,'tag','m_pcolor'));
delete(findobj(ax,'tag','stations'));
delete(findobj(ax,'tag','events'));

% update map
held=ishold(ax);
if(~held); hold(ax,'on'); end
mmap('image',{arf.latlon(:,:,1) arf.latlon(:,:,2) double(arf.beam)},...
    'st',[arf.stla arf.stlo],'ev',arf.latlon0,'parent',ax);
if(~held); hold(ax,'off'); end

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

end
