function [fh]=p3(data,xlimits,ylimits,scale,style,label,fh,sfh) 
%P3   Plot evenly spaced record section
%
%    Description: Plots timeseries and xy records spaced out evenly in the 
%     same plot.  Other record types are ignored.  Optional inputs are the 
%     x/y plot limits ([low high] - y-axis is record number), the amplitude
%     scaling factor (in y-axis units), the style of normalization ('record'
%     or 'global' - scale each record independently or as a group), the 
%     label logical (y-axis ticks have record names for labels - default is
%     0 which is off) and the figure and subplot handles (to put p3 plot in
%     a particular figure and subplot). Default amplitude scaling is 1/10th
%     the number of records.  Default normalization style is 'record'.  
%     Output is the figure handle.
%
%    Usage:  [fh]=p3(data,xlim,ylim,scale,style,labelnames,fh,sfh)
%
%    Notes:  The use of ylim is there but it really isn't useful as you can
%            just select a specific range of data using data(a:b) to do the
%            same thing.
%
%    See also:  p1, p2, recsec

% check number of inputs
error(nargchk(1,7,nargin))

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
set(gcf,'Name','P3 -- SAC Seismogram Plotting Utility', ...
    'NumberTitle','off','color',bgc,'Pointer','crosshair');
set(gca,'FontName',tn,'FontWeight',tw,'FontSize',ts,'box','on',...
    'xcolor',fgc,'ycolor',fgc,'TickDir','out','ticklength',[0 0]);
xlabel('Time (sec)');
ylabel('Record');

% number of records
nrecs=length(data);

% defaults
if(nargin<4 || isempty(scale)); scale=nrecs/10; end
if(nargin<5 || isempty(style)); style='record'; end
if(~any(strcmp(style,{'record' 'global'})))
    error('bad normalization style')
end

% record coloring
colors=hsv(nrecs);

% header info
[b,npts,delta,leven,iftype]=gh(data,'b','npts','delta','leven','iftype');
vers=unique([data.version]);
nver=length(vers);
h(nver)=sachi(vers(nver));
for i=1:nver-1
    h(i)=sachi(vers(i));
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
    plot(time,i+data(i).x*scale/ampmax,'color',colors(i,:),'linewidth',lw);
    plotted(i)=1;
end
hold off

% zooming
axis tight;
if(nargin>1 && ~isempty(xlimits)); axis auto; xlim(xlimits); end
if(nargin>2 && ~isempty(ylimits)); ylim(ylimits); end

% label with record names
if(nargin>5 && label==1)
    if(isfield(data,'name'))
        set(gca,'ytick',1:nrecs);
        set(gca,'yticklabel',{data.name}.');
    end
end

% title (total plotted / total records)
title([num2str(nnz(plotted)) '/' num2str(nrecs) ' Records']);

end