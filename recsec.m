function [fh,lh]=recsec(data,xlimits,ylimits,scale,style,legend_ok,fh,sfh) 
%RECSEC   Plots SAClab data records in a distance spaced record section
%
%    Description: Plots timeseries and xy SAClab records spaced out by the 
%     'gcarc' header field.  Other record types are ignored.  Optional
%     inputs are the x/y plot limits (y-axis is gcarc), the amplitude
%     scaling factor in y-axis units, the style of normalization ('record'
%     or 'global' - scale each record independently or as a group), the
%     legend logical (default is false - no legend), and the figure and
%     subplot handles (to put recsec plot in a particular figure and 
%     subplot).  The legend can be edited and moved using the mouse. 
%     Default amplitude scaling is 1/10th the total 'gcarc' range.  Default
%     normalization style is 'record' (normalizes each independently).  
%     Outputs are the figure and legend handles.
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
set(gcf,'Name',['RECSEC -- ' P.NAME], ...
    'NumberTitle',P.NUMBERTITLE,...
    'color',P.BGCOLOR,...
    'Pointer',P.POINTER);
set(gca,'FontName',P.FONTNAME,'FontWeight',P.FONTWEIGHT,...
    'FontSize',P.FONTSIZE,'box',P.BOX,...
    'xcolor',P.FGCOLOR,'ycolor',P.FGCOLOR,...
    'TickDir',P.TICKDIR,'ticklength',P.TICKLEN);
xlabel('Time (sec)');
ylabel('Distance (\circ)');
grid(P.GRID);

% number of records
nrecs=length(data);

% record coloring
cmap=str2func(P.COLORMAP);
colors=cmap(nrecs);

% header info
leven=glgc(data,'leven');
error(lgcchk('leven',leven))
iftype=genumdesc(data,'iftype');
[b,npts,delta,gcarc]=gh(data,'b','npts','delta','gcarc');

% defaults
if(nargin<4 || isempty(scale))
    scale=(max(gcarc)-min(gcarc))/10;
    
    % fix for no range case
    if(scale==0); scale=1; end 
end
if(nargin<5 || isempty(style)); style='record'; end
if(~any(strcmp(style,{'record' 'global'})))
    warning('SAClab:recsec:badInput','bad normalization style')
    style='record';
end

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
    plot(time,gcarc(i)+data(i).x*scale/ampmax,...
        'color',colors(i,:),'linewidth',P.TRACEWIDTH);
end
hold off

% zooming
axis(P.AXIS);
if(nargin>1 && ~isempty(xlimits)); axis auto; xlim(xlimits); end
if(nargin>2 && ~isempty(ylimits)); ylim(ylimits); end

% legend
if(nargin>5 && legend_ok==1)
    if(isfield(data,'name')) 
        lh=legend(data(goodfiles).name,'location','best');
    else
        lh=legend(strcat({'Record '},...
            cellstr(num2str(indices.')))); 
    end
    set(lh,'interpreter','none')
end

% title (total plotted / total records)
title([num2str(length(indices)) '/' num2str(nrecs) ' Records']);

end
