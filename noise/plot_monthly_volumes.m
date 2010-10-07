function [fh,ax,sx,sy,st,sc]=plot_monthly_volumes(name,cmp,dblim,zerodb)
%PLOT_MONTHLY_VOLUMES    Makes 4x3 grid of monthly fk spectra plots
%
%    Usage:    [fh,ax,sx,sy,st,scb]=plot_monthly_volumes(name,cmp,dblim,zerodb)
%
%    Description:
%     [FH,AX,SX,SY,ST,SCB]=PLOT_MONTHLY_VOLUMES(NAME,CMP,DBLIM,ZERODB)
%
%    Notes:
%
%    Examples:
%
%    See also: MAKE_MONTHLY_HORZ_VOLUMES, MAKE_MONTHLY_Z_VOLUMES

% todo:
% - need to standardize/generalize the inputs
%   - band (should be in Hz) - allow inputting our own, but have a default
%   - remove net, name code below and just use arrayname input ('GAMSEIS')
% - remove polish section from creation
%   - should be a separate function that we pass the figure handle
%     and it just figures out the axes and gets things right
%   - this allows personalized versions of polish (may need variations per matlab version)
%   - makes this operation far less time consuming since we don't have to
%     create the plots every time


% bands
bands=[100 50;
        50 40;
        40 30;
        30 25;
        25 20;
        20 15;
        15 10;
        10 7.5;
        7.5 5];
fbands=1./bands;

% network
switch lower(name)
    case {'gam' 'gamseis'}
        name='Gamseis';
    case {'xb' 'cam' 'cameroon'}
        name='Cameroon';
    case {'xd' 'tan' 'tanzania'}
        name='Tanzania';
    case {'xa' 'sa' 'southafrica'}
        name='South Africa';
    case {'xi' 'eth' 'ethiopia'}
        name='Ethiopia';
    otherwise
        error('bad name');
end

% initialize plots and axes
nb=size(bands,1);
fh=nan(nb,1); ax=cell(nb,1); tax=ax;
for i=1:nb
    fh(i)=figure('color','w','name',...
        [name '   ' num2str(bands(i,2)) 's-' num2str(bands(i,1)) 's']);
    ax{i}=makesubplots(4,3,[],'align','parent',fh(i));
    tax{i}=ax{i}';
    drawnow;
end

% read in monthly volumes and plot
month={'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' ...
    'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};
for i=1:12
    % read in volume
    vol=load(['svol.' lower(cmp) '.' num2str(i,'%02d') '.mat']);
    %vol=fkstructfix(vol);
    
    % loop over each band
    for j=1:nb
        % get fkmap
        map=fkvol2map(vol,fbands(j,:));
        
        % plot it
        plotfkmap(map,dblim,zerodb,'k','w',ax{j}(i));
        
        % fix title
        title(ax{j}(i),month{i});
    end
    drawnow;
end

% save figures
for i=1:nb
    saveas(fh(i),[lower(name) '_' lower(cmp) '_' num2str(bands(i,2)) ...
        's-' num2str(bands(i,1)) 's_' zerodb '_' num2str(dblim(1)) ...
        'db-' num2str(dblim(2)) 'db_temp.fig'],'fig');
end

% polish
nocolorbars([ax{:}]);
nolabels([ax{:}],'xylabel');
[sx,sy,st,sc]=deal(nan(nb,1));
drawnow;
for i=1:nb
    % expand to a portrait
    fig2print(fh(i),'tall');
    drawnow;
    axexpand(ax{i},75);
    drawnow;
    
    % drop extra labels
    nolabels(ax{i}(:,1:3),'xtick');
    nolabels(ax{i}(2:3,:),'ytick');
    drawnow;
    
    % add super labeling
    sx(i)=superxlabel(ax{i},'East/West Slowness (s/deg)',...
        'fontweight','bold','fontsize',12);
    sy(i)=superylabel(ax{i},'North/South Slowness (s/deg)',...
        'fontweight','bold','fontsize',12);
    st(i)=supertitle(ax{i},[name '   ' ...
        num2str(bands(i,2)) 's-' num2str(bands(i,1)) 's'],...
        'fontweight','bold','fontsize',12);
    sc(i)=supercolorbar(ax{i},'location','southoutside');
    set(sc(i),'fontsize',10,'fontweight','bold');
    xlabel(sc(i),['dB (relative to ' zerodb ')'],...
        'fontsize',12,'fontweight','bold');
    drawnow;
    
    % refine position of super labels and plot titles
    set(sx(i),'position',get(sx(i),'position')+[0 0.05 0]);
    drawnow;
    scpos=get(sc(i),'position');
    set(sc(i),'position',...
        [scpos(1) scpos(2)+scpos(4)*0.25 scpos(3) scpos(4)*.5]);
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
        [lower(name) '_' lower(cmp) '_' num2str(bands(i,2)) ...
        's-' num2str(bands(i,1)) 's_' zerodb '_' num2str(dblim(1)) ...
        'db-' num2str(dblim(2)) 'db.pdf']);
    saveas(fh(i),[lower(name) '_' lower(cmp) '_' num2str(bands(i,2)) ...
        's-' num2str(bands(i,1)) 's_' zerodb '_' num2str(dblim(1)) ...
        'db-' num2str(dblim(2)) 'db.fig'],'fig');
end

end
