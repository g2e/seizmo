function [varargout]=adjust_yrmo_plots(fh,mytitle)
%ADJUST_YRMO_PLOTS    Adjust year-month fk spectra figures for printing
%
%    Usage:    adjust_yrmo_plots(fh)
%              adjust_yrmo_plots(fh,mytitle)
%              [ax,sx,sy,st,sc]=adjust_yrmo_plots(...)
%
%    Description:
%     ADJUST_YRMO_PLOTS(FH) adjusts figures by PLOT_YRMO_VOLUMES for
%     printing and saves them to the current directory.  This will probably
%     require tweaking on your part.  FH is a figure handle or array of
%     figure handles.
%
%     ADJUST_YRMO_PLOTS(FH,MYTITLE) uses MYTITLE for the super title for
%     each figure in FH.
%
%     [AX,SX,SY,ST,SC]=ADJUST_YRMO_PLOTS(...) returns the axes and super
%     labeling handles.
%
%    Notes:
%     - Adjusted figures are saved as:
%        fkyrmo_CMP_YEAR_BANDLOs-BANDHIs_ZERODB_DBLIM1db-DBLIM2db.fig
%        fkyrmo_CMP_YEAR_BANDLOs-BANDHIs_ZERODB_DBLIM1db-DBLIM2db.pdf
%
%    Examples:
%     % Create fk volumes, create plots, and prepare for printing:
%     make_yrmo_z_volumes(stack_dir);
%     fh=plot_yrmo_volumes;
%     adjust_yrmo_plots(fh,'My Array Name');
%
%    See also: PLOT_YRMO_VOLUMES, MAKE_YRMO_HORZ_VOLUMES,
%              MAKE_YRMO_Z_VOLUMES

%     Version History:
%        Feb. 15, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 15, 2011 at 16:05 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% default title
if(nargin<2); mytitle=''; end

% check input
if(any(~ishandle(fh(:))) && strcmp(unique(get(fh(:),'tag')),'fkyrmo'))
    error('seizmo:adjust_yrmo_plots:badInput',...
        'FH must be handles to figures made by PLOT_YRMO_VOLUMES!');
end
if(~isstring(mytitle))
    error('seizmo:adjust_yrmo_plots:badInput',...
        'MYTITLE must be a string!');
end

% loop over each figure
nb=numel(fh);
ax=cell(nb,1);
[sx,sy,st,sc]=deal(nan(nb,1));
for i=1:nb
    % raise figure
    figure(fh(i));
    
    % retrieve fkmap plot handles
    % - note we cheated by tagging unplotted axes too
    ax{i}=findobj(fh(i),'tag','fkmap');
    
    % put in proper order
    ax{i}=flipud(ax{i});
    
    % rearrange axes into 3x4
    % - this could be troublesome
    % - we could sort by position to be sure
    ax{i}=reshape(ax{i},3,4);
    
    % remove colorbars
    nocolorbars(ax{i});
    drawnow;
    
    % remove x- & y-axis labels
    nolabels(ax{i},'xylabel');
    drawnow;
    
    % expand figure to a portrait print
    fig2print(fh(i),'tall');
    drawnow;
    
    % drop extra labels
    [tmp,tmp]=max(fliplr(vis),[],2);
    lowermost=midx2li([3 4],5-tmp,2); % right
    [tmp,tmp]=max(vis,[],1);
    leftmost=midx2li([3 4],tmp,1); % top
    nolabels(ax{i}(~ismember(1:9,lowermost)),'xtick');
    nolabels(ax{i}(~ismember(1:9,leftmost)),'ytick');
    drawnow;
    
    % get band
    userdata=get(fh(i),'userdata');
    year=num2str(userdata.year);
    band=userdata.band;
    cmp=userdata.cmp;
    dblim=userdata.dblim;
    zerodb=userdata.zerodb;
    
    % expand plots
    axexpand(ax{i},75);
    drawnow;
    
    % add super labeling
    sx(i)=superxlabel(ax{i},'East/West Slowness (s/deg)',...
        'fontweight','bold','fontsize',12);
    sy(i)=superylabel(ax{i},'North/South Slowness (s/deg)',...
        'fontweight','bold','fontsize',12);
    st(i)=supertitle(ax{i},[mytitle '   ' cmp '   ' year '   ' ...
        num2str(1./band(2)) 's-' num2str(1./band(1)) 's'],...
        'fontweight','bold','fontsize',12);
    set(fh(i),'name',[mytitle '   ' cmp '   ' year '   ' ...
        num2str(1./band(2)) 's-' num2str(1./band(1)) 's']);
    sc(i)=supercolorbar(ax{i},'location','southoutside');
    set(sc(i),'fontsize',10,'fontweight','bold');
    xlabel(sc(i),['dB (relative to ' zerodb ')'],...
        'fontsize',12,'fontweight','bold');
    drawnow;
    
    % refine position of super labels and plot titles
    set(sx(i),'position',get(sx(i),'position')+[0 0.05 0]);
    drawnow;
    scpos=get(sc(i),'position');
    set(sc(i),...
        'position',[scpos(1) scpos(2)+scpos(4)*0.25 scpos(3) scpos(4)*.5]);
    set(sc(i),'xaxislocation','bottom');
    drawnow;
    kids1=get(fh(i),'children');
    titles=get(kids1,'title'); % title handles
    for j=1:numel(titles)
        set(titles{j},'units','normalized');
        set(titles{j},'position',get(titles{j},'position')-[0 0.03 0]);
        drawnow;
    end
    
    % final tweaking
    axstretch(kids1,'v',-15);
    drawnow;
    axmove(kids1,0.05,0);
    drawnow;
    
    % save to pdf & fig
    print(['-f' num2str(fh(i))],'-dpdf',...
        ['fkyrmo_' lower(cmp) '_' year '_' num2str(1/band(2)) ...
        's-' num2str(1/band(1)) 's_' zerodb '_' num2str(dblim(1)) ...
        'db-' num2str(dblim(2)) 'db.pdf']);
    saveas(fh(i),['fkyrmo_' lower(cmp) '_' year '_' num2str(1/band(2)) ...
        's-' num2str(1/band(1)) 's_' zerodb '_' num2str(dblim(1)) ...
        'db-' num2str(dblim(2)) 'db.fig'],'fig');
end

if(nargout)
    varargout={ax,sx,sy,st,sc};
end

end
