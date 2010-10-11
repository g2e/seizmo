function [grp,ax]=selectclusters(data,grp,opt,varargin)
%SELECTCLUSTERS    Select clustered SEIZMO data records graphically
%
%    Usage:    [grp,ax]=selectclusters(data,grp)
%              [grp,ax]=selectclusters(data,grp,option)
%              [grp,ax]=selectclusters(data,grp,option,'option',value,...)
%
%    Description:
%     [GRP,AX]=SELECTCLUSTERS(DATA,GRP) returns the clusters that are
%     graphically selected by the user.  Clusters are defined by the struct
%     GRP (see USERCLUSTER for more details).  The associated waveforms are
%     in SEIZMO struct DATA.  Selection/unselection of clusters is done by
%     left-clicking over a cluster.  Complete the cluster selection by
%     middle-clicking over a cluster or closing the figure.  Outputs are
%     the logical array SELECTED which indicates the clusters that were
%     selected and AX gives the plot handle(s).
%
%     [GRP,AX]=SELECTCLUSTERS(DATA,GRP,OPTION) sets whether selected
%     clusters are kept (good) or deleted (bad).  OPTION must be either
%     'keep' or 'delete'.  When OPTION is 'keep', the background color for
%     selected records is set to a dark green.  For OPTION set to 'delete',
%     the background color is set to a dark red for selected clusters.  The
%     default is 'keep'.
%
%     [GRP,AX]=SELECTCLUSTERS(DATA,GRP,OPTION,'OPTION',VALUE,...) passes
%     plotting options (all arguments after SELECTED) to PLOTCLUSTERS.
%
%    Notes:
%
%    Examples:
%     % Only allow selection of the first three clusters:
%     selectclusters(data,grp,'keep','clusters',1:3)
%
%     % Alter the cluster pre-selection (no kept/good):
%     grp.good=false;
%     grp=selectclusters(data,grp);
%
%    See also: PLOTCLUSTERS, USERCLUSTER, SELECTRECORDS

%     Version History:
%        Sep. 18, 2010 - initial version
%        Sep. 21, 2010 - altered inputs/outputs
%        Oct. 10, 2010 - minor warning fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 10, 2010 at 20:00 GMT

% todo:

% check nargin
if(nargin<1)
    error('seizmo:selectclusters:notEnoughInputs',...
        'Not enough input arguments.');
elseif(nargin>3 && ~mod(nargin,2))
    error('seizmo:selectclusters:plotOptionMustBePaired',...
        'Plot options must be paired with a value!');
end

% check data structure
versioninfo(data,'dep');

% number of records
nrecs=numel(data);

% check cluster struct
if(~isstruct(grp) ...
        || any(~ismember({'T' 'color' 'good'},fieldnames(grp))) ...
        || numel(grp.T)~=nrecs ...
        || ~isequal([max(grp.T) 3],size(grp.color)))
    error('seizmo:selectclusters:badInput',...
        'GRP must be a struct with fields T, good & color!');
end

% number of groups
ngrp=max(grp.T);

% selected
selected=grp.good;

% valid values for options
valid.OPT={'keep' 'delete' 'good' 'bad'};

% default inputs
if(nargin<3 || isempty(opt)); opt='keep'; end

% check inputs
if(~isstring(opt) || isempty(strmatch(lower(opt),valid.OPT)))
    error('seizmo:selectclusters:badInput',...
        ['OPT must be one of the following strings:\n' ...
        sprintf('''%s'' ',valid.OPT{:})]);
elseif((isnumeric(selected) && (any(selected~=fix(selected)) ...
        || any(selected<1 | selected>ngrp))) || (islogical(selected) ...
        && ~any(numel(selected)~=[1 ngrp])))
    error('seizmo:selectclusters:badInput',...
        'GRP.GOOD must be TRUE, FALSE, logical array or linear indices!');
end

% fix selected
if(islogical(selected) && isscalar(selected))
    selected(1:ngrp,1)=selected;
elseif(isnumeric(selected))
    lidx=selected;
    selected=false(ngrp,1);
    selected(lidx)=true;
end

% set color and flip selected if delete
keep=~isempty(strmatch(opt,{'keep' 'good'}));
if(keep)
    color=[0 .3 0];
else % delete
    color=[.3 0 0];
    selected=~selected;
end

% plot clusters
[ax,shown]=plotclusters(data,grp,varargin{:});

% color preselected
bgcolors=get(ax,'color');
if(iscell(bgcolors))
    bgcolors=cell2mat(bgcolors);
end
set(ax(selected(shown)),'color',color);

% loop until user is done
button=0;
while(button~=2)
    % get mouse button pressed
    try
        [x,y,button]=ginput(1);
    catch
        ax=-1;
        button=2;
    end
    if(button==1)
        % grab axis handle
        handle=gca;

        % figure out which record
        clicked=shown(find(handle==ax,1));
        if(isempty(clicked)); continue; end

        % remove from list if in list and change color
        if(selected(clicked))
            selected(clicked)=false;
            set(handle,'color',bgcolors(clicked,:));
            % otherwise add to list and change color
        else
            selected(clicked)=true;
            set(handle,'color',color);
        end
    end
end

% push selected into struct
if(keep)
    grp.good=selected;
else % delete
    grp.good=~selected;
end

end
