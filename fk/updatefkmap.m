function [varargout]=updatefkmap(map,ax)
%UPDATEFKMAP    Quickly updates an existing fk plot with a new map
%
%    Usage:    updatefkmap(map,ax)
%              ax=updatefkmap(map)
%
%    Description: UPDATEFKMAP(MAP,AX) draws an new fk response map given by
%     MAP in an existing axes AX.  This is mainly intended for exploring 3D
%     and 4D fk datasets and for making movies in a faster fashion.
%
%     AX=UPDATEFKMAP(MAP) is the same as calling PLOTFKMAP(MAP) -- ie. a
%     new figure is drawn.
%
%    Notes:
%
%    Examples:
%     Slide through a few fk responses:
%      spolar(1)=fkmap(data,25,201,[1/50 1/45],true);
%      spolar(2)=fkmap(data,25,201,[1/45 1/40],true);
%      spolar(3)=fkmap(data,25,201,[1/40 1/35],true);
%      spolar(4)=fkmap(data,25,201,[1/35 1/30],true);
%      spolar(5)=fkmap(data,25,201,[1/30 1/25],true);
%      spolar(6)=fkmap(data,25,201,[1/25 1/20],true);
%      ax=plotfkmap(spolar(1));
%      for i=2:6
%          pause(1);
%          updatefkmap(spolar(i),ax);
%      end
%
%    See also: PLOTFKMAP, FKMAP, FKVOLUME, FK4D

%     Version History:
%        May  11, 2010 - initial version
%        May  21, 2010 - display period rather than frequency
%        May  26, 2010 - updated for new plotfkmap args (requires passing
%                        info through userdata)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  26, 2010 at 10:20 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check fk struct
error(chkfkstruct(map));

% don't allow array/volume
if(~isscalar(map) || map.volume)
    error('seizmo:updatefkmap:badInput',...
        'MAP must be a scalar fk struct and not a volume!');
end

% just replot if ax isn't an axes handle
if(nargin<2 || ~isscalar(ax) || ~isreal(ax) || ~ishandle(ax) ...
        || ~strcmp('axes',get(ax,'type')))
    ax=plotfkmap(map);
    if(nargout); varargout{1}=ax; end
    return;
end

% plotting function call depends on polar
if(map.polar)
    updatefkpolarmap(map,ax);
else % cartesian
    updatefkcartmap(map,ax);
end
if(nargout); varargout{1}=ax; end

end


function updatefkpolarmap(map,ax)

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

% rescale response
switch zerodb
    case 'min'
        map.response=map.response-min(map.response(:));
        map.normdb=map.normdb+min(map.response(:));
    case 'max'
        map.response=map.response-max(map.response(:));
        map.normdb=map.normdb+max(map.response(:));
    case 'median'
        map.response=map.response-median(map.response(:));
        map.normdb=map.normdb+median(map.response(:));
    case 'abs'
        map.response=map.response+map.normdb;
        map.normdb=0;
end

axes(ax);
h=get(ax,'children');
delete(h);
hold on
nx=numel(map.x);
ny=numel(map.y);
[x,y]=pol2cart(pi/180*map.x(ones(ny,1),:),map.y(:,ones(nx,1)));
pcolor(x,y,map.response);
shading flat;
hold off
set(get(ax,'Title'),'string',...
    {['Number of Stations:  ' num2str(map.nsta)] ...
    ['Begin Time:  ' sprintf('%d.%03d %02d:%02d:%02g',map.butc) ' UTC'] ...
    ['End Time  :  ' sprintf('%d.%03d %02d:%02d:%02g',map.eutc) ' UTC'] ...
    ['Period Range:    ' num2str(1/min(map.z)) ' to ' ...
    num2str(1/max(map.z)) 'Sec'] ...
    ['0 dB = ' num2str(map.normdb) 'dB'] ...
    '' '' ''})

end


function updatefkcartmap(map,ax)

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

% rescale response
switch zerodb
    case 'min'
        map.response=map.response-min(map.response(:));
        map.normdb=map.normdb+min(map.response(:));
    case 'max'
        map.response=map.response-max(map.response(:));
        map.normdb=map.normdb+max(map.response(:));
    case 'median'
        map.response=map.response-median(map.response(:));
        map.normdb=map.normdb+median(map.response(:));
    case 'abs'
        map.response=map.response+map.normdb;
        map.normdb=0;
end

axes(ax);
h=findobj(ax,'type','image');
delete(h);
hold on
imagesc(map.x,map.y,map.response);
childs=get(ax,'children');
set(ax,'children',[childs(2:end); childs(1)]);
hold off
set(get(ax,'Title'),'string',...
    {['Number of Stations:  ' num2str(map.nsta)] ...
    ['Begin Time:  ' sprintf('%d.%03d %02d:%02d:%02g',map.butc) ' UTC'] ...
    ['End Time  :  ' sprintf('%d.%03d %02d:%02d:%02g',map.eutc) ' UTC'] ...
    ['Period Range:    ' num2str(1/min(map.z)) ' to ' ...
    num2str(1/max(map.z)) 'Sec'] ...
    ['0 dB = ' num2str(map.normdb) 'dB']})

end

