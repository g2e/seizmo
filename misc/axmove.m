function axmove(ax,r,u)
%AXMOVE    Moves a set of axes by the specified amount
%
%    Usage:    axmove(ax,right,up)
%
%    Description:
%     AXMOVE(AX,RIGHT,UP) moves the specified axes to the right & up using
%     RIGHT & UP.  Both RIGHT & UP are in figure normalized units.  None of
%     the axes are resized.
%
%    Notes:
%
%    Examples:
%     % move the current axis up by 10% of the figure height
%      axmove(gca,0,.1);
%
%    See also: AXEXPAND, COMPACTAXES

%     Version History:
%        Aug.  5, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug.  5, 2010 at 12:25 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% defaults
if(isempty(ax)); ax=gca; end
if(isempty(r)); r=zeros(size(ax)); end
if(nargin<3 || isempty(u)); u=zeros(size(ax)); end

% check inputs
if(~isreal(ax) || any(~ishandle(ax(:))) ...
        || any(~strcmp('axes',get(ax(:),'type'))))
    error('seizmo:axmove:badInput',...
        'AX must be a valid axes handle!');
elseif(~isreal(r) || ~isreal(u))
    error('seizmo:axmove:badInput',...
        'RIGHT & UP must be real-valued!');
elseif(~isequalsizeorscalar(ax,r,u))
    error('seizmo:axmove:badInput',...
        'All inputs must be equal sized or scalar!');
end

% expand inputs
[ax,r,u]=expandscalars(ax,r,u);

% loop over each axis
for i=1:numel(ax)
    % current position
    pos=get(ax(i),'position');
    
    % new position
    set(ax(i),'position',[pos(1)+r(i) pos(2)+u(i) pos(3:4)]);
end

end
