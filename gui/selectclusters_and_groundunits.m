function [data0,grp,units,ax]=selectclusters_and_groundunits(data,grp,varargin)
%SELECTCLUSTERS_AND_GROUNDUNITS    Work out response issues for clusters
%
%    Usage:    [data,grp,units,ax]=selectclusters_and_groundunits(data,grp)
%              [data,grp,units,ax]=selectclusters_and_groundunits(...
%                         data,grp,...,'option',value,...)
%
%    Description:
%     [DATA,GRP,UNITS,AX]=SELECTCLUSTERS_AND_GROUNDUNITS(DATA,GRP) allows
%     the user to both select clusters and alter the ground units (to check
%     if there are simple response issues -- like records that are actually
%     in velocity rather than in displacement).  To alter the units, the
%     user must right click and then press either "+" or "-" to integrate
%     or differentiate.  Left-clicking a cluster will toggle selection and
%     middle-clicking will exit the gui.  UNITS is a NRECSx1 matrix of
%     integer values ranging from -2 to 2 that indicate how many times the
%     data was/needs to be integrated -- negative numbers indicate
%     differentiation.
%
%     [DATA,GRP,UNITS,AX]=SELECTCLUSTERS_AND_GROUNDUNITS(...
%     DATA,GRP,...,'OPTION',VALUE,..) passes additional options to
%     PLOTCLUSTERS.  See PLOTCLUSTERS for details.
%
%    Notes:
%
%    Examples:
%
%    See also: PLOTCLUSTERS, USERCLUSTER, SELECTRECORDS, SELECTCLUSTERS,
%              ADJUSTCLUSTERS, INTEGRATE, DIFFERENTIATE, REMOVESACPZ

%     Version History:
%        Dec. 17, 2010 - initial version
%        Jan.  6, 2011 - use key2zoompan
%        Jan. 16, 2011 - fix data not matching choice if we alter the units
%                        and then go back to the original (note plot did
%                        not even show correct data)
%        Apr.  3, 2012 - use seizmocheck
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 20:00 GMT

% todo:

% check nargin
if(nargin<1)
    error('seizmo:selectclusters_and_groundunits:notEnoughInputs',...
        'Not enough input arguments.');
elseif(nargin>3 && ~mod(nargin,2))
    error('seizmo:selectclusters_and_groundunits:plotOptionNotPaired',...
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
    error('seizmo:selectclusters_and_groundunits:badInput',...
        'GRP must be a struct with fields T, good & color!');
end

% number of groups
ngrp=max(grp.T);

% selected
selected=grp.good;

% selection setup
keep=true;
color=[0 .3 0];

% check inputs
if((isnumeric(selected) && (any(selected~=fix(selected)) ...
        || any(selected<1 | selected>ngrp))) || (islogical(selected) ...
        && ~any(numel(selected)~=[1 ngrp])))
    error('seizmo:selectclusters_and_groundunits:badInput',...
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
gunits=zeros(ngrp,1);
data0=data; % expensive ...
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
    clicked=shown(find(handle==ax,1));
    if(isempty(clicked)); continue; end
    
    % act based on button
    if(button==1) % toggle selection
        % remove from list if in list and change color
        if(selected(clicked))
            selected(clicked)=false;
            set(handle,'color',bgcolors(clicked,:));
        % otherwise add to list and change color
        else
            selected(clicked)=true;
            set(handle,'color',color);
        end
    elseif(button==43) % integrate
        % skip if already double int
        if(gunits(clicked)==2); continue; end
        
        % integrate cluster data and replot them
        gunits(clicked)=gunits(clicked)+1;
        recs=grp.T==clicked;
        switch gunits(clicked)
            case -2 % double diff
                data0(recs)=differentiate(differentiate(data(recs)));
            case -1 % diff
                data0(recs)=differentiate(data(recs));
            case 0 % nada
                data0(recs)=data(recs);
            case 1 % int
                data0(recs)=integrate(data(recs));
            case 2 % double int
                data0(recs)=integrate(integrate(data(recs)));
        end
        plotclusters(data0,grp,varargin{:},...
            'clusters',clicked,'axis',handle);
    elseif(button==45) % differentiate
        % skip if already double diff
        if(gunits(clicked)==-2); continue; end
        
        % differentiate cluster data and replot them
        gunits(clicked)=gunits(clicked)-1;
        recs=grp.T==clicked;
        switch gunits(clicked)
            case -2 % double diff
                data0(recs)=differentiate(differentiate(data(recs)));
            case -1 % diff
                data0(recs)=differentiate(data(recs));
            case 0 % nada
                data0(recs)=data(recs);
            case 1 % int
                data0(recs)=integrate(data(recs));
            case 2 % double int
                data0(recs)=integrate(integrate(data(recs)));
        end
        plotclusters(data0,grp,varargin{:},...
            'clusters',clicked,'axis',handle);
    else
        key2zoompan(button,handle);
    end
end

% group units to individual units
units=zeros(size(grp.T));
for i=1:max(grp.T)
    units(grp.T==i)=gunits(i);
end

% push selected into struct
if(keep)
    grp.good=selected;
else % delete
    grp.good=~selected;
end

end
