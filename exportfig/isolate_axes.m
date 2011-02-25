%ISOLATE_AXES Isolate the specified axes in a figure on their own
%
% Examples:
%   fh = isolate_axes(ah)
%   fh = isolate_axes(ah, vis)
%
% This function will create a new figure containing the axes specified, and
% also their associated legends and colorbars. The axes specified must all
% be in the same figure, but they will generally only be a subset of the
% axes in the figure.
%
% IN:
%    ah - An array of axes handles, which must come from the same figure.
%    vis - A boolean indicating whether the new figure should be visible.
%          Default: false.
%
% OUT:
%    fh - The handle of the created figure.

% Copyright (C) Oliver Woodford 2011

function fh = isolate_axes(ah, vis)
% Make sure we have an array of handles
if ~all(ishandle(ah))
    error('ah must be an array of handles');
end
% Check that the handles are all for axes, and are all in the same figure
fh = ancestor(ah(1), 'figure');
nAx = numel(ah);
for a = 1:nAx
    if ~strcmp(get(ah(a), 'Type'), 'axes')
        error('All handles must be axes handles.');
    end
    if ~isequal(ancestor(ah(a), 'figure'), fh)
        error('Axes must all come from the same figure.');
    end
end
% Tag the axes so we can find them in the copy
old_tag = get(ah, 'Tag');
if nAx == 1
    old_tag = {old_tag};
end
set(ah, 'Tag', 'ObjectToCopy');
% Create a new figure exactly the same as the old one
fh = copyfig(fh); %copyobj(fh, 0);
if nargin < 2 || ~vis
    set(fh, 'Visible', 'off');
end
% Reset the axes tags
for a = 1:nAx
    set(ah(a), 'Tag', old_tag{a});
end
% Get all the axes and gui objects
axs = get(fh, 'Children');
% Find the objects to save
ah = findobj(axs, 'Tag', 'ObjectToCopy');
if numel(ah) ~= nAx
    error('Incorrect number of axes found.');
end
I = ~ismember(axs, ah);
% Set the axes tags to what they should be
for a = 1:nAx
    set(ah(a), 'Tag', old_tag{a});
end
% Keep any legends and colorbars which overlap the subplots
lh = findobj(axs, 'Type', 'axes', '-and', {'Tag', 'legend', '-or', 'Tag', 'Colorbar'})';
nLeg = numel(lh);
if nLeg > 0
    [leg_ind leg_ind] = ismember(lh, axs); 
    ax_pos = get(ah, 'OuterPosition');
    if nAx > 1
        ax_pos = cell2mat(ax_pos(:));
    end
    ax_pos(:,3:4) = ax_pos(:,3:4) + ax_pos(:,1:2);
    leg_pos = get(lh, 'OuterPosition');
    if nLeg > 1;
        leg_pos = cell2mat(leg_pos(:));
    end
    leg_pos(:,3:4) = leg_pos(:,3:4) + leg_pos(:,1:2);
    for a = 1:nAx
            % Overlap test
            I(leg_ind(leg_pos(:,1) < ax_pos(a,3) & leg_pos(:,2) < ax_pos(a,4) &...
                      leg_pos(:,3) > ax_pos(a,1) & leg_pos(:,4) > ax_pos(a,2))) = false;
    end
end
% Delete all axes except for the input axes and associated items
delete(axs(I));
return

function fh = copyfig(fh)
% Is there a legend?
if isempty(findobj(fh, 'Type', 'axes', 'Tag', 'legend'))
    % Safe to copy using copyobj
    fh = copyobj(fh, 0);
else
    % copyobj will change the figure, so save and then load it instead
    tmp_nam = [tempname '.fig'];
    hgsave(fh, tmp_nam);
    fh = hgload(tmp_nam);
    delete(tmp_nam);
end
return