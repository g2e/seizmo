function [fh,sfh]=p1(data,xlimits,ylimits,ncols,fh) 
%P1    Plot SAClab data records in individual subplots
%
%    Description: Plots timeseries and xy records in individual subplots.
%     Other filetypes are ignored (a space is left in the figure though).
%     Flags for 'o','a','f','t(n)' are drawn if defined.  Labels are drawn
%     for the flags using 'ko','ka','kf','kt(n)' if defined or the field
%     name (o,a,f,t(n)).  Optional inputs are limits for x and y axis, the
%     number of columns to use in the figure, and the figure handle. Output
%     is the figure and subfigure handles.
% 
%    Usage:  [fh,sfh]=p1(data,xlim,ylim,ncols,fh)
%
%    See also:  p2, recsec

% check number of inputs
error(nargchk(1,5,nargin))

% check data structure
if(~isstruct(data))
    error('input data is not a structure')
elseif(~isvector(data))
    error('data structure not a vector')
elseif(~isfield(data,'version') || ~isfield(data,'head') || ...
        ~isfield(data,'x'))
    error('data structure does not have proper fields')
end

% default basic plotting parameters
bgc='k';    % background color
fgc='w';    % axis color
ts=4;       % text size
tn='times'; % text name
tw='light'; % text weight
lw=1;       % line width of record
fw=3;       % flag width

% initialize plot
if(nargin<5 || isempty(fh) || fh<1); fh=figure;
else figure(fh); end
whitebg(bgc);
set(gcf,'Name','P1 -- SAC Seismogram Plotting Utility', ...
    'NumberTitle','off','color',bgc,'Pointer','crosshair');

% header info
[t,kt,o,ko,a,ka,f,kf,b,delta,npts,gcarc,leven,iftype]=...
    gh(data,'t','kt','o','ko','a','ka','f','kf','b','delta','npts',...
    'gcarc','leven','iftype');
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end

% trim strings
kt=strtrim(kt);
ko=strtrim(ko);
ka=strtrim(ka);
kf=strtrim(kf);

% number of records
nrecs=length(data);

% record coloring
colors=hsv(nrecs);

% default columns/rows
if(nargin<4 || isempty(ncols)); ncols=fix(sqrt(nrecs)); nrows=ceil(nrecs/ncols);
else nrows=ceil(nrecs/ncols); end

% check user defined x/y limits (force into individual subplot control)
if(nargin>1 && ~isempty(xlimits))
    if(isvector(xlimits))
        xlimits=xlimits(:).';
        if(length(xlimits)~=2); xlimits=[];
        else xlimits=xlimits(ones(nrecs,1),:);
        end
    elseif(any(size(xlimits)==nrecs))
        if(size(xlimits,2)==nrecs); xlimits=xlimits.'; end
        if(all(size(xlimits)~=2)); xlimits=[]; end
    else xlimits=[];
    end
end
if(nargin>2 && ~isempty(ylimits))
    if(isvector(ylimits))
        ylimits=ylimits(:).';
        if(length(ylimits)~=2); ylimits=[];
        else ylimits=ylimits(ones(nrecs,1),:);
        end
    elseif(any(size(ylimits)==nrecs))
        if(size(ylimits,2)==nrecs); ylimits=ylimits.'; end
        if(all(size(ylimits)~=2)); ylimits=[]; end
    else ylimits=[];
    end
end

% loop through each record
sfh=zeros(nrecs,1);
for i=1:nrecs
    % header version
    v=data(i).version==vers;
    
    % check if timeseries or xy (if not skip)
    if(~any(iftype(i)==[h(v).enum(1).val.itime h(v).enum(1).val.ixy])); continue; end
    
    % get record timing
    if(leven(i)==h(v).true); time=b(i)+(0:npts(i)-1).'*delta(i);
    else time=data(i).t; end
    
    % focus to new subplot and draw record
    sfh(i)=subplot(nrows,ncols,i);
    plot(time,data(i).x,'color',colors(i,:),'linewidth',lw);
    set(gca,'FontName',tn,'FontWeight',tw,'FontSize',ts, ...
        'xcolor',fgc,'ycolor',fgc, 'TickDir','out','ticklength',[0 0]);
    grid on;
    
    % zooming
    axis tight;
    if(nargin>1 && ~isempty(xlimits)); axis auto; xlim(xlimits(i,:)); end
    if(nargin>2 && ~isempty(ylimits)); ylim(ylimits(i,:)); end
    
    % scaling parameters for header flags
    fylimits=get(gca,'ylim');
    fxlimits=get(gca,'xlim');
    yrange=fylimits(2)-fylimits(1);
    xrange=fxlimits(2)-fxlimits(1);
    ypad=0.05*yrange;  % these are simple - really you will have to
    xpad=0.02*xrange;  % adjust these based on your font and plot size
    
    % plot origin flag (orange)
    hold on
    if(o(i)~=h(v).undef.ntype)
        plot([o(i) o(i)].',[fylimits(1)+ypad fylimits(2)-2*ypad].',...
            'color',[1 0.5 0],'linewidth',fw);
        if(~strcmp(ko{i},h(v).undef.stype))
            text(o(i)+xpad,fylimits(2)-2*ypad,ko{i},'color',[1 0.5 0], ...
                'fontname',tn,'fontweight',tw,'verticalalignment','top',...
                'clipping','on','fontsize',ts);
        else
            text(o(i)+xpad,fylimits(2)-2*ypad,'o','color',[1 0.5 0], ...
                'fontname',tn,'fontweight',tw,'verticalalignment','top',...
                'clipping','on','fontsize',ts);
        end
    end
    
    % plot arrival flag (green)
    if(a(i)~=h(v).undef.ntype)
        plot([a(i) a(i)].',[fylimits(1)+ypad fylimits(2)-2*ypad].',...
            'color','g','linewidth',fw);
        if(~strcmp(ka{i},h(v).undef.stype))
            text(a(i)+xpad,fylimits(2)-2*ypad,ka{i},'color','g', ...
                'fontname',tn,'fontweight',tw,'verticalalignment','top',...
                'clipping','on','fontsize',ts);
        else
            text(a(i)+xpad,fylimits(2)-2*ypad,'a','color','g', ...
                'fontname',tn,'fontweight',tw,'verticalalignment','top',...
                'clipping','on','fontsize',ts);
        end
    end
    
    % plot finish flag (red)
    if(f(i)~=h(v).undef.ntype)
        plot([f(i) f(i)].',[fylimits(1)+ypad fylimits(2)-2*ypad].',...
            'color','r','linewidth',fw);
        if(~strcmp(kf{i},h(v).undef.stype))
            text(f(i)+xpad,fylimits(2)-2*ypad,kf{i},'color','r', ...
                'fontname',tn,'fontweight',tw,'verticalalignment','top',...
                'clipping','on','fontsize',ts);
        else
            text(f(i)+xpad,fylimits(2)-2*ypad,'f','color','r', ...
                'fontname',tn,'fontweight',tw,'verticalalignment','top',...
                'clipping','on','fontsize',ts);
        end
    end
    
    % plot picks (yellow)
    for j=0:9
        if(t(i,j+1)~=h(v).undef.ntype)
            plot([t(i,j+1) t(i,j+1)].',[fylimits(1)+ypad fylimits(2)-(1+mod(j,5))*ypad].',...
                'color','y','linewidth',fw)
            if(~strcmp(kt{i,j+1},h(v).undef.stype))
                text(t(i,j+1)+xpad,fylimits(2)-(1+mod(j,5))*ypad,kt{i,j+1},'color','y', ...
                    'fontname',tn,'fontweight',tw,'fontsize',ts,...
                    'verticalalignment','top','clipping','on');
            else
                text(t(i,j+1)+xpad,fylimits(2)-(1+mod(j,5))*ypad,['t' num2str(j)],...
                    'color','y','fontname',tn,'fontweight',tw,...
                    'verticalalignment','top','clipping','on',...
                    'fontsize',ts);
            end
        end
    end
    
    % record name and degree distance
    if(isfield(data,'name') && ~isempty(data(i).name))
        title([texlabel(data(i).name,'literal') '  -  ' num2str(gcarc(i)) '\circ'])
    else title(['RECORD ' num2str(i) '  -  ' num2str(gcarc(i)) '\circ']); 
    end
    hold off
end

end
