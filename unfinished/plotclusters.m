function [h,sh]=plotclusters(T,pop,big,small,amps,h)
%PLOTCLUSTERS    Plots interactive population cutoff and grouped traces subplots

% GLOBALS
global CONF GCOLOR

% NUMBER OF SAMPLES/GROUPS
ns=size(amps,1);
ng=max(T);
c1=length(big);
c2=length(small);

% PLOT SMALL GROUPS
if(c2)
    % NEW PLOT
    if(nargin==5)
        h(2)=figure('numbertitle','off','name','SMALL GROUPS',...
            'menubar','none','color',CONF.BGCOLOR,'Pointer','crosshair');
    end
    
    % SMALL GROUPS SUBPLOT SETUP
    ncol=ceil(sqrt(c2));
    nrow=ceil(c2/ncol);
    
    % FIGURE SETUP
    figure(h(2));
    whitebg(CONF.BGCOLOR)
    set(gcf,'numbertitle','off','name','SMALL GROUPS',...
        'menubar','none','color',CONF.BGCOLOR,'Pointer','crosshair');
    
    % GROUP SUBPLOTS
    for i=1:c2
        subplot(nrow,ncol,i)
        plot((0:ns-1)*CONF.DELTA+CONF.NOWWIN(1),amps(:,small(i)==T),...
            'color',GCOLOR(small(i),:))
        set(gca,'fontweight',CONF.FONTWEIGHT,'fontsize',CONF.FONTSIZE,...
            'FontName',CONF.FONTNAME,'xcolor',CONF.FGCOLOR,...
            'ycolor',CONF.FGCOLOR)
        hold on
        text(CONF.NOWWIN(1)+ns*CONF.DELTA*0.05,0.5,...
            ['GROUP ' num2str(small(i))],...
            'fontsize',CONF.FONTSIZE,'fontweight',CONF.FONTWEIGHT)
        hold off
        axis tight
        grid on
    end
    drawnow;
end

% PLOT BIG GROUPS NEW PLOT
if(nargin==5)
    h(1)=figure('numbertitle','off','name','BIG GROUPS',...
        'menubar','none','color',CONF.BGCOLOR,'Pointer','crosshair');
end

% BIG GROUPS SUBPLOT SETUP
ncol=max([1 ceil(sqrt(c1))]);
nrow=ceil(c1/ncol)+1;

% FIGURE SETUP
figure(h(1));
whitebg(CONF.BGCOLOR)
set(gcf,'numbertitle','off','name','BIG GROUPS',...
    'menubar','none','color',CONF.BGCOLOR,'Pointer','crosshair');

% GROUP POPULATION PLOT
sh=subplot(nrow,ncol,1:ncol);
plot([1:ng; 1:ng],[pop.'; ones(1,ng)],'-k','linewidth',2)
hold on
scatter(1:ng,pop,100,GCOLOR,'h','filled','MarkerEdgeColor','k')
plot([0 max(T)+1],[CONF.MINSIG CONF.MINSIG],'r','linewidth',4)
set(gca,'fontweight',CONF.FONTWEIGHT,'fontsize',CONF.FONTSIZE,...
    'FontName',CONF.FONTNAME,'xcolor',CONF.FGCOLOR,'ycolor',CONF.FGCOLOR,...
    'yscale','log','xtick',big,'ygrid','on',...
    'xaxislocation','top','ticklength',[0 0])
text(0,CONF.MINSIG,'GROUP POP CUT',...
    'fontweight',CONF.FONTWEIGHT,'fontsize',CONF.FONTSIZE,...
    'verticalalignment','bottom')
if(CONF.FILTER)
    title(CONF.TITLE,...
        'fontweight',CONF.FONTWEIGHT,'fontsize',CONF.FONTSIZE)
else
    title('UNFILTERED','fontweight',CONF.FONTWEIGHT,...
        'fontsize',CONF.FONTSIZE)
end
ylabel('POPULATION','fontweight',CONF.FONTWEIGHT,'fontsize',CONF.FONTSIZE)
xlim([0 max(T)+1])
hold off

% GROUP SUBPLOTS
for i=1:c1
    subplot(nrow,ncol,ncol+i)
    plot((0:ns-1)*CONF.DELTA+CONF.NOWWIN(1),amps(:,big(i)==T),...
        'color',GCOLOR(big(i),:))
    set(gca,'fontweight',CONF.FONTWEIGHT,'fontsize',CONF.FONTSIZE,...
        'FontName',CONF.FONTNAME,'xcolor',CONF.FGCOLOR,'ycolor',CONF.FGCOLOR)
    hold on
    text(CONF.NOWWIN(1)+ns*CONF.DELTA*0.05,0.5,...
        ['GROUP ' num2str(big(i))],...
        'fontsize',CONF.FONTSIZE,'fontweight',CONF.FONTWEIGHT)
    hold off
    axis tight
    grid on
end

end
