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
%     GRP (see USERCLUSTER for struct details).  The associated waveforms
%     are in SEIZMO struct DATA.  Selection/unselection of clusters is
%     performed by left-clicking a cluster.  Complete the cluster selection
%     by middle-clicking over a cluster or closing the figure.  Outputs are
%     the modified struct GRP which indicates the clusters that were
%     selected in field GRP.good and AX gives the plot handle(s).
%
%     [GRP,AX]=SELECTCLUSTERS(DATA,GRP,OPTION) sets whether selected
%     clusters are kept (good) or deleted (bad).  OPTION must be either
%     'keep' or 'delete'.  When OPTION is 'keep', the background color for
%     selected records is set to a dark green.  For OPTION set to 'delete',
%     the background color is set to a dark red for selected clusters.  The
%     default is 'keep'.  Output GRP.good indicates the "kept" clusters
%     with a corresponding TRUE element while "deleted" clusters are FALSE.
%
%     [GRP,AX]=SELECTCLUSTERS(DATA,GRP,OPTION,'PARAMETER',VALUE,...) passes
%     plotting options (all arguments after OPTION) to PLOTCLUSTERS.
%
%    Notes:
%     - Preselection is done using the input GRP.good values
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
%        Jan.  6, 2011 - use key2zoompan
%        Jan. 13, 2011 - doc cleanup, handle plotcluster no plot case, fix
%                        bugs in cases with empty clusters
%        Mar. 31, 2011 - fix bug in checking of grp.good
%        Apr.  3, 2012 - use seizmocheck
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 20:00 GMT

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
error(seizmocheck(data,'dep'));

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
        && ~any(numel(selected)==[1 ngrp])))
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

% check if no plot
if(isempty(ax))
    % push pre-selected into struct
    if(keep)
        grp.good=selected;
    else % delete
        grp.good=~selected;
    end
    return;
end

% color preselected
bgcolors=get(ax,'color');
if(iscell(bgcolors))
    bgcolors=cell2mat(bgcolors);
end
set(ax(ismember(shown,find(selected))),'color',color);

% loop until user is done
button=0;
while(button~=2)
    % get mouse button pressed
    try
        [x,y,button]=ginput(1);
    catch
        ax=-1;
        break;
    end
    
    % grab axis handle
    handle=gca;
    
    % figure out which record
    clicked=shown(handle==ax);
    if(isempty(clicked)); continue; end
    
    % act based on button
    if(button==1)
        % remove from list if in list and change color
        if(selected(clicked))
            selected(clicked)=false;
            set(handle,'color',bgcolors(handle==ax,:));
        % otherwise add to list and change color
        else
            selected(clicked)=true;
            set(handle,'color',color);
        end
    else
        key2zoompan(button,handle);
    end
end

% push selected into struct
if(keep)
    grp.good=selected;
else % delete
    grp.good=~selected;
end

end
