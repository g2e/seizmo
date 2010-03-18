function [fh]=population(grp,fh,sh)
%POPULATION    Makes a stem plot of cluster populations

CONF.FGCOLOR='k';
CONF.FONTNAME='arial';
CONF.FONTWEIGHT='bold';
CONF.FONTSIZE=12;
CONF.MINSIG=3;

ng=max(grp.T);
pop=histc(grp.T,1:ng);
big=find(pop>=CONF.MINSIG);
GCOLOR=grp.color;


%figure(fh);
%subplot(sh);
figure;
set(gca,'colororder',GCOLOR,'nextplot','replacechildren')
plot([1:ng; 1:ng],[pop.'; 0.9*ones(1,ng)],'-.','linewidth',3.5)
hold on
scatter(1:ng,pop,200,GCOLOR,'p','filled','MarkerEdgeColor','k')
plot([0 ng+1],[CONF.MINSIG CONF.MINSIG],'r','linewidth',4)
set(gca,'fontweight',CONF.FONTWEIGHT,'fontsize',CONF.FONTSIZE,...
    'FontName',CONF.FONTNAME,'xcolor',CONF.FGCOLOR,'ycolor',CONF.FGCOLOR,...
    'yscale','log','xtick',big,'ygrid','on','ticklength',[0 0])
text(0,CONF.MINSIG,'GROUP POP CUT',...
    'fontweight',CONF.FONTWEIGHT,'fontsize',CONF.FONTSIZE,...
    'verticalalignment','bottom')
ylabel('POPULATION','fontweight',CONF.FONTWEIGHT,'fontsize',CONF.FONTSIZE)
xlim([0 ng+1])
ylim([0.9 max(pop)*1.1])
hold off

end
