function [varargout]=movekids(kids,move)
%MOVEKIDS     Moves the specified children plot objects
%
%    Usage:    movekids(kids)
%              movekids(kids,position)
%              kids=movekids(...)
%
%    Description:
%     MOVEKIDS(KIDS) places the child objects specified by the handles in
%     KIDS to the TOP/FRONT of their parent object visually.  All handles
%     in KIDS must share the same parent.
%
%     MOVEKIDS(KIDS,POSITION) alters where the child objects are placed.
%     Valid position strings are 'top'/'front' or 'back'/'bottom'.  The
%     default is 'top' or 'front' (they do the same action), which places
%     the specified children objects at the front.
%
%     KIDS=MOVEKIDS(...) returns the reordered child object handles for the
%     parent of the input KIDS.
%
%    Notes:
%     - MOVEKIDS is a simplified version of UISTACK.
%
%    Examples:
%     % bring an object in the back of the current axes to the front
%     kids=get(gca,'children');
%     kids=movekids(kids(end),'top');
%
%    See also: UISTACK

%     Version History:
%        Aug. 26, 2010 - initial version
%        Sep. 14, 2010 - handle arrays of handles better
%        Nov. 20, 2010 - better siblings check
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 20, 2010 at 20:50 GMTs

% todo:

% check nargin
error(nargchk(1,2,nargin));

% defaults (no move)
if(nargin<2 || isempty(move)); move='top'; end
%if(nargin<3 || isempty(dist)); dist=0; end

% check kids (are handles, same parent)
kids=kids(:); % force as column vector
if(~isreal(kids) || any(~ishandle(kids)))
     error('seizmo:movekids:kidsOver12',...
         'Some of your kids are not!');
end

% check kids are siblings
if(isscalar(kids))
    p=get(kids,'parent');
else
    p=unique(cell2mat(get(kids,'parent')));
    if(~isscalar(p))
        error('seizmo:movekids:notSiblings',...
            'KIDS must be siblings!');
    end
end

% check move/dist
if(~ischar(move) || ndims(move)~=2 || size(move,1)~=1)
    error('seizmo:movekids:badMove',...
        'MOVE must be a string!');
end
%if(~isreal(dist) || ~isscalar(dist) || dist~=fix(dist))
%    error('seizmo:movekids:badDistance',...
%        'DIST must be a scalar integer!');
%end

% parent and other kids
allkids=get(p,'children');
bastards=~ismember(allkids,kids);
%nk=numel(allkids);

% act by move type
idx=strmatch(lower(move),{'top' 'bottom' 'front' 'back'});
if(isempty(idx))
    error('seizmo:movekids:badMove',...
        'Unknown move: %s !',move);
end
switch idx
    case {1 3} % top/front
        set(p,'children',[kids; allkids(bastards)]);
    case {2 4} % bottom/back
        set(p,'children',[allkids(bastards); kids]);
    otherwise
        error('seizmo:movekids:badMove',...
            'Unknown move: %s !',move);
end

% output kids
if(nargout); varargout{1}=get(p,'children'); end

end
