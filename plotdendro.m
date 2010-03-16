function [perm,colors,fh,sfh]=plotdendro(data,Z,varargin)
%PLOTDENDRO    Plots correlation linkage and seismograms
%
%    Usage: [PERM,COLORS,FH,SFH]=...
%             plotdendro(Z,DATA,'plot_option',plot_option_value,...)
%
%    Description: Plots the heirarchial linkage given in Z corresponding to
%     records given in DATA.  Records in DATA are plotted with respect to
%     their placement in the tree for enhanced analysis.  Clustering and
%     cluster colors can be controlled with 'treelimit', 'treecolormap',
%     'treelimitcolor' and 'treedefcolor'.  To switch the xaxis to
%     similarity rather than dissimilarity set 'dissim' to false.  Outputs
%     are the permutation vector PERM that gives ordering of records in the
%     tree.  COLORS are the individual colors of the records permuted with
%     PERM.  FH and SFH give the figure and subplot handles.
%
%    Examples:
%     Get similarities using mcxc, assemble the tree using linkage, and
%     visualize the heirarchy using pdendro:
%       records=records2mat(data);
%       cv=mcxc(records,'npeaks',1);
%       Z=linkage(1-cv.','average');
%       plotdendro(Z,data,'treelimit',1);
%
%    See also: PLOT1, PLOT2, PLOT0, RECORDSECTION, MCXC,
%              DENDROGRAM, LINKAGE, CLUSTER (Statistics Toolbox)

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% get plotting style defaults
P=plotconfig;

% allow access to plot styling using global SEIZMO structure
global SEIZMO; fields=fieldnames(P).';
for i=fields; if(isfield(SEIZMO,i)); P.(i{:})=SEIZMO.(i{:}); end; end

% sort out optional arguments
for i=1:2:length(varargin)
    % plotting parameter field
    varargin{i}=upper(varargin{i});
    if(isfield(P,varargin{i}))
        P.(varargin{i})=varargin{i+1};
    else
        warning('seizmo:pdendro:badInput',...
            'Unknown Option: %s',varargin{i}); 
    end
end

% clean up unset parameters
P=plotconfigfix(P);

% number of records
nrecs=numel(data);

% select/open plot
if(isempty(P.FIGHANDLE) || P.FIGHANDLE<1)
    fh=figure;
else
    fh=figure(P.FIGHANDLE);
end

% style the figure (early)
set(fh,'name',['PDENDRO -- ' P.NAME],...
        'numbertitle',P.NUMBERTITLE,...
        'menubar',P.MENUBAR,...
        'toolbar',P.TOOLBAR,...
        'renderer',P.RENDERER,...
        'renderermode',P.RENDERERMODE,...
        'doublebuffer',P.DOUBLEBUFFER,...
        'color',P.FIGBGCOLOR,...
        'pointer',P.POINTER);

% dendrogram subplot (on the left, descending to the right)
sfh(1)=subplot('Position',[0.1 0.1 0.4 0.8]);
S=cell(1,nrecs); S(:)={''};
try
    % my modified dendrogram (for coloring)
    [H,H2,perm]=ddendrogram(Z,0,'labels',S,...
        'orientation','left','colorthreshold',P.TREELIMIT,...
        'colormap',P.TREECOLORMAP,'defaultcolor',P.TREEDEFCOLOR);
catch
    % Matlab's dendrogram
    [H,H2,perm]=dendrogram(Z,0,'labels',S,...
        'orientation','left','colorthreshold',P.TREELIMIT);
end
perm=perm.';

% style the plot
set(sfh(1),'fontname',P.AXISFONT,'fontweight',P.AXISFONTWEIGHT,...
    'fontsize',P.AXISFONTSIZE,'box',P.BOX,...
    'xscale',P.XSCALE,'yscale',P.YSCALE,...
    'xcolor',P.XAXISCOLOR,'ycolor',P.YAXISCOLOR,...
    'xaxislocation',P.XAXISLOCATION,...
    'tickdir',P.TICKDIR,'ticklength',P.TICKLEN,...
    'linewidth',P.AXISLINEWIDTH,'color',P.PLOTBGCOLOR,...
    'ytick',[],'xminortick','on');
grid(P.GRID);

% dissimilarity
if(P.DISSIM)
    dendrolabel='DISSIMILARITY';
% similarity
else
    set(sfh(1),'xticklabel',1-get(sfh(1),'xtick'))
    set(sfh(1),'xtickmode','manual')
    dendrolabel='SIMILARITY';
end

% label the plot
if(isempty(P.TREETITLE)); P.TREETITLE='LINKAGE DENDOGRAM'; end
if(isempty(P.TREEXLABEL)); P.TREEXLABEL=dendrolabel; end
title(P.TREETITLE,'fontname',P.TITLEFONT,'fontweight',P.TITLEFONTWEIGHT,...
    'fontsize',P.TITLEFONTSIZE,'color',P.TITLEFONTCOLOR,...
    'interpreter',P.TITLEINTERP);
xlabel(P.TREEXLABEL,'fontname',P.XLABELFONT,...
    'fontweight',P.XLABELFONTWEIGHT,...
    'fontsize',P.XLABELFONTSIZE,'color',P.XLABELFONTCOLOR,...
    'interpreter',P.XLABELINTERP);

% retrieve colors
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
clear temp

% waveform subplot
sfh(2)=subplot('Position',[0.5 0.1 0.4 0.8]);
set(sfh(2),'color',P.PLOTBGCOLOR);
plot0(data(perm),'fighandle',fh,'subhandle',sfh(2),...
    'yaxislocation','right',varargin{:},...
    'colormap',colors,'markers',false);

% match dendrogram yaxis to waveform plot
span=get(sfh(2),'YLim');
set(sfh(1),'Ylim',span);

% restyle the figure (late)
set(fh,'name',['PDENDRO -- ' P.NAME],...
        'numbertitle',P.NUMBERTITLE,...
        'menubar',P.MENUBAR,...
        'toolbar',P.TOOLBAR,...
        'renderer',P.RENDERER,...
        'doublebuffer',P.DOUBLEBUFFER,...
        'color',P.FIGBGCOLOR,...
        'pointer',P.POINTER);

% add dissimilarity limit marker
subplot(sfh(1));
hold on
plot([P.TREELIMIT P.TREELIMIT],span,P.TREELIMITCOLOR,...
    'linewidth',P.MARKERWIDTH);
if(P.MARKERLABELS)
    text(P.TREELIMIT,0,'DISSIMILARITY LIMIT','fontname',P.MARKERFONT,...
        'fontweight',P.MARKERFONTWEIGHT,'fontsize',P.MARKERFONTSIZE,...
        'verticalalignment','bottom','rotation',90,...
        'color',P.TREELIMITCOLOR,'clipping','on');
end
hold off

end
