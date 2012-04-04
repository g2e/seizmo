function [varargout]=plotpop(grp,ax)
%PLOTPOP    Makes a stem plot of cluster populations
%
%    Usage:    plotpop(grp)
%              plotpop(grp,ax)
%              ax=plotpop(...)
%
%    Description:
%     PLOTPOP(GRP) makes a stem plot of the cluster populations given in
%     the struct GRP.  See USERCLUSTER for the format of GRP.  This is
%     mainly for interactive population-based cluster cutting using POPCUT.
%
%     PLOTPOP(GRP,AX) sets the axis handle AX to draw in.  This is useful
%     for subplotting.
%
%     AX=PLOTPOP(...) returns the plot handle.
%
%    Notes:
%
%    Examples:
%     % Cluster some data interactively then look at the populations:
%     grp=usercluster(data,cg);
%     plotpop(grp)
%
%    See also: POPCUT, PLOTCLUSTERS, USERCLUSTER

%     Version History:
%        Mar. 23, 2010 - initial version
%        Sep. 18, 2010 - major update
%        Oct.  6, 2010 - handle 0 pop clusters
%        Mar. 31, 2011 - drawnow call avoids buggy matlab scatter issues,
%                        improve ticking
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 17:35 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% defaults
if(nargin<2); ax=[]; end

% check inputs
if(~isstruct(grp) || any(~ismember({'T' 'color'},fieldnames(grp))))
    error('seizmo:plotpop:badInput',...
        'GRP must be a struct with fields T & color!');
end
if(~isempty(ax) && (~isscalar(ax) || ~isreal(ax)))
    error('seizmo:plotpop:badInput',...
        'AX must be a handle to a single axis!');
end

% check the axes
if(isempty(ax) || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')))
    % new figure
    fh=figure('color','k');
    ax=axes('parent',fh);
else
    cla(ax,'reset');
    axes(ax);
end

% cluster info
ngrp=max(grp.T);
pop=histc(grp.T,1:ngrp);
ok=find(pop>0);
ok=ok(:)';
nok=numel(ok);

% style plot
set(ax,'colororder',grp.color,'nextplot','replacechildren');

% plot stems & stars
plot(ax,[ok; ok],[pop(ok).'; 0.9*ones(1,nok)],'-.','linewidth',2);
hold(ax,'on');
scatter(ax,ok,pop(ok),200,grp.color,'p','filled');
drawnow; % avoids matlab scatter bugs
xlabel(ax,'Cluster ID','fontweight','bold','fontsize',10);
ylabel(ax,'Population','fontweight','bold','fontsize',10);
set(ax,'xlim',[0 ngrp+1])
set(ax,'ylim',[0.9 max(pop)*1.1])
hold(ax,'off');

% style plot
set(ax,'yscale','log','ticklength',[0 0],'ygrid','on','box','on',...
    'color','k','xcolor','w','ycolor','w',...
    'fontsize',10,'fontweight','bold');

% manual yticks b/c matlab log ticks are a failure
ylimits=ylim(ax);
high=max(fix(log10(ylimits)));
if(high<1); high=1; set(ax,'ylim',[ylimits(1) 10]); end
set(ax,'ytick',10.^(0:high),'yticklabel',10.^(0:high));

% output
if(nargout); varargout={ax}; end

end
