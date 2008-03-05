function [fh]=p3(data,xlimits,ylimits,scale,style,label,fh,sfh) 
%P3   Plot SAClab data records in an evenly spaced record section
%
%    Description: Plots timeseries and xy SAClab records spaced out evenly
%     in the same plot.  Other record types are ignored.  Optional inputs
%     are the x/y plot limits ([low high] - y-axis is record number), the
%     amplitude scaling factor (in y-axis units), the style of normalizing
%     ('record' or 'global' - scale each record separately or as a group),
%     the label logical (label y-axis ticks with record names - default is
%     false) and the figure and subplot handles (to put p3 plot in a
%     particular figure and subplot). Default amplitude scaling is 1/10th
%     the number of records.  Default normalization style is 'record'.  
%     Output is the figure handle.
%
%    Usage:  [fh]=p3(data,xlim,ylim,scale,style,labelnames,fh,sfh)
%
%    Notes:  The use of ylim is there but it really isn't useful as you can
%            just select a specific range of data using data(a:b) to do the
%            same thing.
%
%    Examples:
%
%    See also:  p1, p2, recsec

% check number of inputs
error(nargchk(1,8,nargin))

% check data structure
error(seischk(data,'x'))

% plotting style defaults
P=pconf;

% allow access to plot styling using global SACLAB structure
global SACLAB
fields=fieldnames(P).';
for i=fields; if(isfield(SACLAB,i)); P.(i{:})=SACLAB.(i{:}); end; end

% initialize plot
if(nargin<7 || isempty(fh) || fh<1); fh=figure;
else figure(fh); if(nargin==8 && ~isempty(sfh)); subplot(sfh); end; end
whitebg(P.BGCOLOR);
set(gcf,'Name',['P3 -- ' P.NAME], ...
    'NumberTitle',P.NUMBERTITLE,...
    'color',P.BGCOLOR,...
    'Pointer',P.POINTER);
set(gca,'FontName',P.FONTNAME,'FontWeight',P.FONTWEIGHT,...
    'FontSize',P.FONTSIZE,'box',P.BOX,...
    'xcolor',P.FGCOLOR,'ycolor',P.FGCOLOR,...
    'TickDir',P.TICKDIR,'ticklength',P.TICKLEN);
grid(P.GRID);

% number of records
nrecs=length(data);

% check options
if(nargin<4 || isempty(scale)); scale=nrecs/10; end
if(nargin<5 || isempty(style)); style='record'; end
if(~any(strcmp(style,{'record' 'global'})))
    warning('SAClab:p3:badInput','bad normalization style')
    style='record';
end

% record coloring
cmap=str2func(P.COLORMAP);
colors=cmap(nrecs);

% header info
leven=glgc(data,'leven');
error(lgcchk('leven',leven))
iftype=genumdesc(data,'iftype');
[b,npts,delta]=gh(data,'b','npts','delta');

% check filetype
goodfiles=(strcmp(iftype,'Time Series File') | ...
    strcmp(iftype,'General X vs Y file'));
indices=find(goodfiles).';

% find max amplitude for global style normalization
if(strcmp(style,'global'));
    ampmax=0;
    for i=indices
        tmax=max(abs(data(i).x));
        if(tmax>ampmax); ampmax=tmax; end
    end
end

% loop through each file
hold on
for i=indices
    % get record timing
    if(strcmp(leven(i),'true')); time=(b(i)+(0:npts(i)-1)*delta(i)).';
    else time=data(i).t; end
    
    % get max amplitude for record style normalization
    if(strcmp(style,'record')); ampmax=max(abs(data(i).x)); end
    
    % plot series
    plot(time,i+data(i).x*scale/ampmax,...
        'color',colors(i,:),'linewidth',P.TRACEWIDTH);
end
hold off

% zooming
axis(P.AXIS);
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
title([num2str(length(indices)) '/' num2str(nrecs) ' Records']);

end
