function [fh,lh]=recsec(data,xlimits,ylimits,scale,style,legend_ok,fh,sfh) 
%RECSEC   Plot distance sorted record section
%
%    Description: Plots timeseries and xy records spaced out by their 
%     'gcarc' header field.  Other record types are ignored.  Optional
%     inputs are the x/y plot limits (y-axis is gcarc), the amplitude
%     scaling factor in y-axis units, the style of normalization ('record'
%     or 'global' - scale each record independently or as a group), the
%     legend logical (default is 0 - no legend), and the figure and subplot
%     handles (to put recsec plot in a particular figure and subplot).  The
%     legend can be edited and moved using the mouse. Default amplitude 
%     scaling is 1/10th the total distance range.  Default normalization 
%     style is 'record'.  Outputs are the figure and legend handles.
%
%    Usage:  [fh,lh]=recsec(data,xlim,ylim,scale,style,legend,fh,sfh)
%
%    Examples:
%     To make a record section of data between 100 and 150 degrees showing 
%     the first 300 seconds with the max amplitude normalized to 3 degrees
%     while preserving the relative amplitudes between records and
%     including a legend
%      recsec(data,[0 300],[100 150],3,'global',1)
%
%    See also:  p1, p2, p3

% check number of inputs
error(nargchk(1,8,nargin))

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

% initialize plot
if(nargin<7 || isempty(fh) || fh<1); fh=figure;
else figure(fh); if(nargin==8 && ~isempty(sfh)); subplot(sfh); end; end
whitebg(bgc);
set(gcf,'Name','RECSEC -- SAC Seismogram Plotting Utility', ...
    'NumberTitle','off','color',bgc,'Pointer','crosshair');
set(gca,'FontName',tn,'FontWeight',tw,'FontSize',ts,'box','on',...
    'xcolor',fgc,'ycolor',fgc,'TickDir','out','ticklength',[0 0]);
xlabel('Time (sec)');
ylabel('Distance (\circ)');
grid on;

% number of records
nrecs=length(data);

% record coloring
colors=hsv(nrecs);

% header info
[b,npts,delta,leven,iftype,gcarc]=gh(data,'b','npts','delta','leven','iftype','gcarc');
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
end

% defaults
if(nargin<4 || isempty(scale)); scale=(max(gcarc)-min(gcarc))/10; end
if(nargin<5 || isempty(style)); style='record'; end
if(~any(strcmp(style,{'record' 'global'})))
    error('bad normalization style')
end

% find max amplitude for global style normalization
if(strcmp(style,'global'));
    ampmax=0;
    for i=1:nrecs
        v=data(i).version==vers;
        if(~any(iftype(i)==[h(v).enum(1).val.itime h(v).enum(1).val.ixy])); continue; end
        tmax=max(abs(data(i).x));
        if(tmax>ampmax); ampmax=tmax; end
    end
end

% loop through each file
hold on
plotted=false(nrecs,1);
for i=1:nrecs
    % header version
    v=data(i).version==vers;
    
    % check if timeseries or xy (if not skip)
    if(~any(iftype(i)==[h(v).enum(1).val.itime h(v).enum(1).val.ixy])); continue; end
    
    % get record timing
    if(leven(i)==h(v).true); time=b(i)+(0:npts(i)-1).'*delta(i);
    else time=data(i).t; end
    
    % get max amplitude for record style normalization
    if(strcmp(style,'record')); ampmax=max(abs(data(i).x)); end
    
    % plot series
    plot(time,gcarc(i)+data(i).x*scale/ampmax,'color',colors(i,:),'linewidth',lw);
    plotted(i)=1;
end
hold off

% zooming
axis tight;
if(nargin>1 && ~isempty(xlimits)); axis auto; xlim(xlimits); end
if(nargin>2 && ~isempty(ylimits)); ylim(ylimits); end

% legend
if(nargin>5 && legend_ok==1)
    if(isfield(data,'name')); lh=legend(data(plotted).name,'location','best');
    else i=1:nrecs; lh=legend(strcat({'Record '},cellstr(num2str(i(plotted).')))); end
    set(lh,'interpreter','none')
    set(lh,'fontsize',4)
end

% title (total plotted / total records)
title([num2str(nnz(plotted)) '/' num2str(nrecs) ' Records']);

end