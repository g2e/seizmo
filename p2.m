function [fh,lh]=p2(data,xlimits,ylimits,legend_ok,norm,fh,sfh) 
%P2    Overlay plot of SAClab data records
%
%    Description: Plots timeseries and xy SAClab records over one another
%     in a plot.  Other record types are ignored.  Optional inputs are the
%     x/y limits for the plot ([low high]), a legend logical (default is
%     false), a normalization flag (default is false), and the figure and 
%     subplot handles (to put p2 in a particular figure and subplot).  The 
%     legend box can be moved by left-clicking on its box and dragging.  
%     The legend labels can be graphically edited by double-clicking on 
%     them.  Outputs are the figure and legend handles.
%
%    Usage:  [fh,lh]=p2(data,xlim,ylim,legend_ok,norm,fh,sfh)
%
%    Examples:
%     To overlay the first 4 records
%      p2(data(1:4)) 
%
%     To overlay the 5th and 8th records from 0 to 300 seconds
%      p2(data([5 8]),[0 300]) 
%
%     To plot all traces with a legend, without limiting the x/y axis
%      p2(data,[],[],1)
%
%    See also:  p1, p3, recsec

% check number of inputs
error(nargchk(1,7,nargin))

% check data structure
error(seischk(data,'x'))

% plotting style defaults
P=pconf;

% allow access to plot styling using global SACLAB structure
global SACLAB
fields=fieldnames(P).';
for i=fields; if(isfield(SACLAB,i)); P.(i{:})=SACLAB.(i{:}); end; end

% initialize plot
if(nargin<6 || isempty(fh) || fh<1); fh=figure;
else figure(fh); if(nargin==7 && ~isempty(sfh)); subplot(sfh); end; end
whitebg(P.BGCOLOR);
set(gcf,'Name',['P2 -- ' P.NAME], ...
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

% record coloring
cmap=str2func(P.COLORMAP);
colors=cmap(nrecs);

% header info
iftype=genumdesc(data,'iftype');
leven=glgc(data,'leven');
error(lgcchk('leven',leven))
[b,npts,delta]=gh(data,'b','npts','delta');

% normalization
scaling=ones(nrecs,1);
if(nargin>4 && ~isempty(norm) && norm); scaling=gnrm(data); end

% loop through each file
hold on
plotted=false(nrecs,1);
for i=1:nrecs
    % check if timeseries or xy (if not skip)
    if(~any(strcmp(iftype(i),...
            {'Time Series File' 'General X vs Y file'})))
        continue; 
    end
    
    % get record timing
    if(strcmp(leven(i),'true')); time=b(i)+(0:npts(i)-1).'*delta(i);
    else time=data(i).t; end
    
    % plot series
    plot(time,data(i).x/scaling(i),...
        'color',colors(i,:),'linewidth',P.TRACEWIDTH);
    plotted(i)=true;
end
hold off

% zooming
axis(P.AXIS);
if(nargin>1 && ~isempty(xlimits)); axis auto; xlim(xlimits); end
if(nargin>2 && ~isempty(ylimits)); ylim(ylimits); end

% legend
if(nargin>3 && ~isempty(legend_ok) && legend_ok)
    if(isfield(data,'name'))
        % filenames
        lh=legend(data(plotted).name,'location','best');
    else
        % record indices
        i=1:nrecs; 
        lh=legend(strcat({'Record '},cellstr(num2str(i(plotted).')))); 
    end
    set(lh,'interpreter','none')
end

% title (total plotted / total records)
title([num2str(nnz(plotted)) '/' num2str(nrecs) ' Records']);

end
