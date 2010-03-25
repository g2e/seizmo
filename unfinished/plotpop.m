function [varargout]=plotpop(grp,h)
%PLOTPOP    Makes a stem plot of cluster populations
%
%    Usage:    plotpop(grp)
%              plotpop(grp,h)
%              h=plotpop(...)
%
%    Description: PLOTPOP(GRP) makes a stem plot of the cluster populations
%     given in the struct GRP.  See USERCLUSTER for the format of GRP.
%     This is mainly for interactive population-based cluster cutting using
%     POPCUT.
%
%     PLOTPOP(GRP,H) sets the figure/plot handle that PLOTPOP should draw
%     in.  Useful for putting this in a subplot.
%
%     H=PLOTPOP(...) returns the figure handle that the plot is in.
%
%    Notes:
%
%    Examples:
%     Cluster some data interactively then look at the populations:
%      grp=usercluster(data,cg);
%      plotpop(grp)
%
%    See also: POPCUT, POPULATION, USERCLUSTER

%     Version History:
%        Mar. 23, 2010 - complete rewrite, renamed from PLOT0 to PLOTEVEN
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 23, 2010 at 17:35 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% check inputs
if(~isstruct(grp) || any(~ismember({'T' 'color'},fieldnames(grp))))
    error('seizmo:plotpop:badInput',...
        'GRP must be a struct with fields T & color!');
end

% plot styling
BGCOLOR='k';
FGCOLOR='w';
FONTNAME='arial';
FONTWEIGHT='bold';
FONTSIZE=12;

% cluster info
NG=max(grp.T);
POP=histc(grp.T,1:NG);
GCOLOR=grp.color;

% pull up figure
if(nargin==2 && isscalar(h) && ishandle(h))
    if(h==fix(h))
        figure(h);
    else
        subplot(h);
        h=get(h,'parent');
    end
else
    h=figure;
end
varargout{1}=h;

% style plot
set(h,'color',BGCOLOR,'name','SEIZMO -- PLOTPOP');
set(gca,'fontname',FONTNAME,'fontweight',FONTWEIGHT,...
    'fontsize',FONTSIZE,'box','on','xcolor',FGCOLOR,'ycolor',FGCOLOR,...
    'colororder',GCOLOR,'nextplot','replacechildren','yscale','log',...
    'ticklength',[0 0],'color',BGCOLOR,'ygrid','on');

% plot stems & stars
plot([1:NG; 1:NG],[POP.'; 0.9*ones(1,NG)],'-.','linewidth',3.5)
hold on
scatter(1:NG,POP,200,GCOLOR,'p','filled','MarkerEdgeColor',FGCOLOR)
xlabel('CLUSTER NUMBER','fontname',FONTNAME,'fontweight',FONTWEIGHT,...
    'fontsize',FONTSIZE,'color',FGCOLOR);
ylabel('POPULATION','fontname',FONTNAME,'fontweight',FONTWEIGHT,...
    'fontsize',FONTSIZE,'color',FGCOLOR);
xlim([0 NG+1])
ylim([0.9 max(POP)*1.1])
hold off

end
