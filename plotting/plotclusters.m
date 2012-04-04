function [varargout]=plotclusters(data,grp,varargin)
%PLOTCLUSTERS    Plots SEIZMO data records as clusters in subplots
%
%    Usage:    plotclusters(data,grp)
%              plotclusters(data,grp,'option',value,...)
%              [ax,clusters]=plotclusters(...)
%
%    Description:
%     PLOTCLUSTERS(DATA,GRP) draws all non-xyz records in SEIZMO struct
%     DATA in separate subplots according to their cluster id in the struct
%     GRP.  GRP should be structured as from USERCLUSTER.
%
%     PLOTCLUSTERS(DATA,GRP,'OPTION',VALUE,...) sets certain plotting
%     options to do simple manipulation of the plots.  Available options
%     are (beyond those available to PLOT2):
%      POPRNG   -- only plot clusters with a population within the given
%                  range (default is [1 inf] - cannot plot 0 pop clusters)
%      CLUSTERS -- only plot clusters in this group
%
%     [AX,CLUSTERS]=PLOTCLUSTERS(...) returns the handles for all the axes
%     drawn in.  This is useful for more detailed plot manipulation.  The
%     2nd output lists the clusters plotted.
%
%    Notes:
%
%    Examples:
%     % From the first 10 clusters only show those
%     % with a population greater than 3:
%     plotclusters(data,grp,'poprng',[3 inf],'clusters',1:10);
%
%    See also: POPCUT, PLOTPOP, USERCLUSTER

%     Version History:
%        Sep. 18, 2010 - initial version
%        Oct.  6, 2010 - handle 0 pop clusters (do not plot)
%        Jan. 13, 2011 - fix no plot varargout bug
%        Apr.  3, 2012 - use seizmocheck
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 17:35 GMT

% todo:

% check nargin
error(nargchk(2,inf,nargin));

% check data struct
error(seizmocheck(data,'dep'));

% number of records
nrecs=numel(data);

% check cluster struct
if(~isstruct(grp) || any(~ismember({'T' 'color'},fieldnames(grp))) ...
        || numel(grp.T)~=nrecs || ~isequal([max(grp.T) 3],size(grp.color)))
    error('seizmo:plotpop:badInput',...
        'GRP must be a struct with fields T & color!');
end

% default/parse options
opt=parse_seizmo_plot_options(varargin{:});

% get populations
ngrp=max(grp.T);
pop=histc(grp.T,1:ngrp);

% default selected clusters
if(isempty(opt.CLUSTERS)); opt.CLUSTERS=1:ngrp; end

% groups to plot
good=find(pop>=opt.POPRNG(1) & pop<=opt.POPRNG(2) & pop>=1);
good=intersect(good,opt.CLUSTERS);
ngrp=numel(good);

% return if none
if(~ngrp); if(nargout); varargout={[] []}; end; return; end

% check axes
if(numel(opt.AXIS)<ngrp || ~isreal(opt.AXIS(1:ngrp)) ...
        || any(~ishandle(opt.AXIS(1:ngrp))) ...
        || any(~strcmp('axes',get(opt.AXIS(1:ngrp),'type'))))
    % new figure
    fh=figure('color',opt.BGCOLOR);
    
    % columns and rows
    if(isempty(opt.NUMCOLS))
        opt.NUMCOLS=fix(sqrt(ngrp));
    end
    nrows=ceil(ngrp/opt.NUMCOLS);
    
    % subplots
    opt.AXIS=makesubplots(nrows,opt.NUMCOLS,1:ngrp,...
        opt.ALIGN{:},'parent',fh);
end

% loop over clusters
th=nan(ngrp,1);
for i=1:ngrp
    % plot cluster
    opt.AXIS(i)=plot2(data(grp.T==good(i)),varargin{:},...
        'cmap',grp.color(good(i),:),'axis',opt.AXIS(i));
    
    % add text
    x=get(opt.AXIS(i),'xlim');
    y=get(opt.AXIS(i),'ylim');
    th(i)=text(x(2)-.01*(x(2)-x(1)),y(2)-.01*(y(2)-y(1)),...
        num2str(good(i)),'color',grp.color(good(i),:),...
        'parent',opt.AXIS(i),...
        'fontsize',2*opt.FONTSIZE,'fontweight',opt.FONTWEIGHT,...
        'horizontalalignment','right','verticalalignment','top');
end

% output
if(nargout); varargout={opt.AXIS(1:ngrp) good}; end

end
