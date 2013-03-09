function [perm,colors,ax]=plotdendro(data,Z,varargin)
%PLOTDENDRO    Plots correlation linkage and seismograms
%
%    Usage: [perm,colors,ax]=plotdendro(data,z)
%           [...]=plotdendro(data,z,'option',value,...)
%
%    Description:
%     [PERM,COLORS,AX]=PLOTDENDRO(DATA,Z) plots the heirarchial linkage
%     of dissimilarities in Z corresponding to records given in SEIZMO
%     struct DATA.  Records in DATA are plotted with respect to their
%     placement in the tree for easy visual inspection.  The clustering
%     cutoff is set to a dissimilarity of 0.2.  PERM gives the ordering of
%     the records in the dendrogram & waveform plots.  COLORS are the
%     individual colors of the permuted records.  AX contains the axes
%     handles to the dendrogram and waveform plots.
%
%     [...]=PLOTDENDRO(DATA,Z,'OPTION',VALUE,...) allows changing certain
%     plotting options to do simple manipulation of the plots as well as
%     controlling the clustering cutoff.  Available options are:
%      FGCOLOR      -- foreground color (axes, text, labels)
%      BGCOLOR      -- background color (does not set figure color)
%      AXIS         -- axes to plot in (need 2)
%      COLORMAP     -- colormap for coloring data
%      XLABEL       -- x axis label
%      YLABEL       -- y axis label
%      TITLE        -- title
%      XLIM         -- x axis limits (tight by default)
%      YLIM         -- y axis limits (tight by default)
%      LINEWIDTH    -- line width of records (default is 1)
%      LINESTYLE    -- line style of records (can be char/cellstr array)
%      NUMCOLS      -- number of subplot columns
%      UTC          -- plot in absolute time if TRUE (UTC, no leap support)
%      DATEFORMAT   -- date format used if ABSOLUTE (default is auto)
%      NORMSTYLE    -- normalize 'individually' or as a 'group'
%      NORMMAX      -- max value of normalized records
%      NORM2YAXIS   -- scale to yaxis range (NORMMAX is fraction of range)
%      NAMESONYAXIS -- true/false or 'kstnm' 'stcmp' 'kname'
%      XDIR         -- 'normal' or 'reverse'
%      YDIR         -- 'normal' or 'reverse'
%      FONTSIZE     -- size of fonts in the axes
%      XSCALE       -- 'linear' or 'log'
%      YSCALE       -- 'linear' or 'log'
%      AMPSCALE     -- 'linear' or 'log'
%      DISTCUT      -- 0 to 1 (default of 0.2)
%      DISTCUTCOLOR -- default is 'r' (red)
%      OTHERCOLOR   -- color of unclustered records/nodes ([.5 .5 .5])
%
%    Notes:
%
%    Examples:
%     % Get similarities using correlate, assemble the tree using linkage,
%     % and visualize the heirarchy using plotdendro:
%     peaks=correlate(data,'mcxc','noauto','normxc','peaks');
%     Z=linkage(1-peaks.cg.','average');
%     plotdendro(data,Z,'cutoff',1);
%
%    See also: PLOT0, PLOT1, PLOT2, RECORDSECTION, CORRELATE,
%              DENDROGRAM, LINKAGE, CLUSTER (Statistics Toolbox)

%     Version History:
%        Aug. 15, 2010 - rewrite
%        Aug. 24, 2010 - require two axes handles
%        Aug. 26, 2010 - cutoff/cutoffcolor options (ie renamed them)
%        Sep. 21, 2010 - distcut/distcutcolor
%        Oct.  1, 2010 - better axis handling
%        Apr. 19, 2011 - fixed cutoff text for plot0 ydir change
%        Jan. 30, 2013 - update example for new correlate
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 30, 2013 at 23:00 GMT

% todo:

% check nargin
error(nargchk(2,inf,nargin));

% check struct
error(seizmocheck(data,'dep'));
nrecs=numel(data);

% check Z
if(~isreal(Z) || ~isequal(size(Z),[nrecs-1 3]))
    error('seizmo:plotdendro:badInput',...
        'Linkage matrix Z corrupt or does not correspond to DATA!');
end

% default/parse options
opt=parse_seizmo_plot_options(varargin{:});

% all in one plot
if(isempty(opt.AXIS) || numel(opt.AXIS)~=2 || ~isreal(opt.AXIS) ...
        || any(~ishandle(opt.AXIS)) ...
        || any(~strcmp('axes',get(opt.AXIS,'type'))))
    % new figure
    figure('color',opt.BGCOLOR);
    ax(2)=subplot('position',[0.5 0.1 0.4 0.8]);
    ax(1)=subplot('position',[0.1 0.1 0.4 0.8]);
else
    % take axis, get pos, delete, make two inside
    ax=opt.AXIS;
    axes(ax(1));
end

% dendrogram subplot (note it uses gca)
S=cell(1,nrecs); S(:)={''}; % no labels
try
    % my modified dendrogram (for coloring)
    [H,H2,perm]=ddendrogram(Z,0,'labels',S,...
        'orientation','left','colorthreshold',opt.CUTOFF,...
        'colormap',opt.CMAP,'defaultcolor',opt.OTHERCOLOR);
catch
    % Matlab's dendrogram
    [H,H2,perm]=dendrogram(Z,0,'labels',S,...
        'orientation','left','colorthreshold',opt.CUTOFF);
end
perm=perm.';

% polish dendrogram plot
set(H,'linewidth',opt.LINEWIDTH);
set(ax(1),'xcolor',opt.FGCOLOR,'ycolor',opt.FGCOLOR,'color',opt.BGCOLOR,...
    'fontsize',opt.FONTSIZE,'ytick',[],'xminortick','on',...
    'ydir',opt.YDIR,'linewidth',opt.LINEWIDTH);
box(ax(1),'on');
grid(ax(1),'on');

% label dendrogram plot
title(ax(1),'LINKAGE DENDOGRAM',...
    'color',opt.FGCOLOR,'fontsize',opt.FONTSIZE);
xlabel(ax(1),'DISSIMILARITY',...
    'color',opt.FGCOLOR,'fontsize',opt.FONTSIZE);

% now retrieve colors from the dendrogram in a ridiculous way
colors=zeros(nrecs,3);
for i=nrecs-1:-1:1
    temp=get(H(i),'ydata');
    if (temp(1)==round(temp(1)))
        colors(temp(1),:)=get(H(i),'color');
    end
    if (temp(4)==round(temp(4)))
        colors(temp(4),:)=get(H(i),'color');
    end
end

% waveform subplot (note waveforms are permuted)
plot0(data(perm),varargin{:},'ax',ax(2),'colormap',colors);
set(ax(2),'yaxislocation','right');

% sync y axes
span=get(ax(2),'ylim');
linkaxes(ax([2 1]),'y');

% add cluster cutoff bar
hold(ax(1),'on');
plot(ax(1),opt.CUTOFF*ones(1,2),span,...
    opt.CUTOFFCOLOR,'linewidth',opt.LINEWIDTH);
text(opt.CUTOFF,span(1)+0.02*diff(span),'DISSIMILARITY LIMIT',...
    'fontsize',opt.FONTSIZE,'color',opt.CUTOFFCOLOR,...
    'verticalalignment','bottom','horizontalalignment','left',...
    'rotation',90,'parent',ax(1));
hold(ax(1),'off');

% put dendrogram on top so interactivity isn't impaired
ph=get(ax(1),'parent');
kids=get(ph,'children');
idx1=find(kids==ax(1));
idx2=find(kids==ax(2));
kids([idx1 idx2])=kids([idx2 idx1]);
set(ph,'children',kids);

end

