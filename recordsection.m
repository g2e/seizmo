function [fh,lh]=recordsection(data,varargin) 
%RECORDSECTION   Plots SEIZMO data records in a distance spaced record section
%
%    Description: Plots timeseries and xy SEIZMO records spaced out by the 
%     'gcarc' header field.  Other record types are ignored.  Optional
%     inputs should correspond to fields returned by function pconf.
%     Outputs are the figure and legend handles.
%
%    Usage:  [fh,lh]=recordsection(data,'plot_option',plot_option_value,...)
%
%    Examples:
%     To make a record section of data between 100 and 150 degrees showing 
%     the first 300 seconds with the max amplitude normalized to 3 degrees
%     while preserving the relative amplitudes between records and
%     including a legend:
%      recsec(data,'xlimits',[0 300],'ylimits',[100 150],...
%           'normstyle','group','norm2yrange',false','normmax',3,...
%           'legend',true)
%
%     To plot records against azimuth/backazimuth:
%       recordsection(data,'yfield','az')
%       recordsection(data,'yfield','baz')
%
%    See also:  plot1, plot2, plot0

% check data structure
error(seizmocheck(data,'dep'))

% get plotting style defaults
P=plotconfig;

% allow access to plot styling using global SEIZMO structure
global SEIZMO; fields=fieldnames(P).';
for i=fields; if(isfield(SEIZMO,i)); P.(i{:})=SEIZMO.(i{:}); end; end

% sort out optional arguments
for i=1:2:length(varargin)
    % plotting parameter field
    varargin{i}=upper(varargin{i});
    if(isfield(P,varargin{i}))
        P.(varargin{i})=varargin{i+1};
    else
        warning('seizmo:recsec:badInput','Unknown Option: %s',varargin{i}); 
    end
end

% clean up unset parameters
P=plotconfigfix(P);

% select/open plot
if(isempty(P.FIGHANDLE) || P.FIGHANDLE<1)
    fh=figure;
else
    fh=figure(P.FIGHANDLE);
    if(~isempty(P.SUBHANDLE)) 
        subplot(P.SUBHANDLE); 
    end; 
end

% SOME STYLING OF THE PLOT
set(gcf,'name',['RECSEC -- ' P.NAME],...
        'numbertitle',P.NUMBERTITLE,...
        'menubar',P.MENUBAR,...
        'toolbar',P.TOOLBAR,...
        'renderer',P.RENDERER,...
        'renderermode',P.RENDERERMODE,...
        'doublebuffer',P.DOUBLEBUFFER,...
        'color',P.FIGBGCOLOR,...
        'pointer',P.POINTER);
set(gca,'fontname',P.AXISFONT,'fontweight',P.AXISFONTWEIGHT,...
        'fontsize',P.AXISFONTSIZE,'box',P.BOX,...
        'xcolor',P.XAXISCOLOR,'ycolor',P.YAXISCOLOR,...
        'xaxislocation',P.XAXISLOCATION,'yaxislocation',P.YAXISLOCATION,...
        'tickdir',P.TICKDIR,'ticklength',P.TICKLEN,...
        'linewidth',P.AXISLINEWIDTH,'color',P.PLOTBGCOLOR,...
        'xminortick',P.XMINORTICK,'yminortick',P.YMINORTICK);
grid(P.GRID);

% number of records
nrecs=length(data);

% record coloring
try
    % try using as a colormap function
    cmap=str2func(P.COLORMAP);
    colors=cmap(nrecs);
catch
    % otherwise its a color array/string
    colors=repmat(P.COLORMAP,ceil(nrecs/size(P.COLORMAP,1)),1);
end

% header info
leven=getlgc(data,'leven');
iftype=getenumdesc(data,'iftype');
[b,npts,delta,depmin,depmax,yfield]=...
    getheader(data,'b','npts','delta','depmin','depmax',P.YFIELD);

% yaxis scaling for amplitudes
if(P.NORM2YRANGE)
    scale=(max(yfield)-min(yfield))*P.NORMMAX;
    
    % fix for no range case
    if(scale==0); scale=1; end
else
    scale=P.NORMMAX;
end

% check normalization style
if(~any(strcmpi(P.NORMSTYLE,{'single' 'group'})))
    warning('seizmo:recsec:badInput','bad normalization style')
    P.NORMSTYLE='single';
end

% check filetype (only timeseries or xy)
goodfiles=(strcmp(iftype,'Time Series File') | ...
    strcmp(iftype,'General X vs Y file'));
indices=find(goodfiles).';

% find amplitude scaling
ampmax=max(abs([depmin.'; depmax.';]));
if(strcmpi(P.NORMSTYLE,'group')); ampmax(:)=max(ampmax(indices)); end

% loop through each file
hold on
for i=indices
    % get record timing
    if(strcmp(leven(i),'true')); time=(b(i)+(0:npts(i)-1)*delta(i)).';
    else time=data(i).ind; end
    
    % plot series
    plot(time,yfield(i)+data(i).dep/ampmax(i)*scale,...
        'color',colors(i,:),'linewidth',P.RECWIDTH);
end
hold off

% zooming
axis(P.AXIS{:});
if(~isempty(P.XLIMITS)); axis auto; xlim(P.XLIMITS); end
if(~isempty(P.YLIMITS)); ylim(P.YLIMITS); end

% legend
if(P.LEGEND)
    if(isfield(data,'name')); lh=legend(data(goodfiles).name);
    else lh=legend(strcat({'Record '},cellstr(num2str(indices.')))); end
    set(lh,'location',P.LEGENDLOCATION,...
        'edgecolor',P.LEGENDBOXCOLOR,...
        'linewidth',P.LEGENDBOXWIDTH,...
        'textcolor',P.LEGENDFONTCOLOR,...
        'box',P.LEGENDBOX,...
        'color',P.LEGENDBGCOLOR,...
        'interpreter',P.LEGENDINTERP,...
        'fontname',P.LEGENDFONT,...
        'fontsize',P.LEGENDFONTSIZE,...
        'fontweight',P.LEGENDFONTWEIGHT);
end

% built in default axis labels/title
% title (total plotted / total records)
if(isempty(P.TITLE))
    P.TITLE=[num2str(length(indices)) '/' num2str(nrecs) ' Records']; 
end
if(isempty(P.XLABEL)); P.XLABEL='Time (sec)'; end
if(isempty(P.YLABEL) && strcmpi(P.YLABEL,'gcarc'))
    P.YLABEL='Distance (\circ)'; 
end
title(P.TITLE,'fontname',P.TITLEFONT,'fontweight',P.TITLEFONTWEIGHT,...
    'fontsize',P.TITLEFONTSIZE,'color',P.TITLEFONTCOLOR,...
    'interpreter',P.TITLEINTERP);
xlabel(P.XLABEL,'fontname',P.XLABELFONT,'fontweight',P.XLABELFONTWEIGHT,...
    'fontsize',P.XLABELFONTSIZE,'color',P.XLABELFONTCOLOR,...
    'interpreter',P.XLABELINTERP);
ylabel(P.YLABEL,'fontname',P.YLABELFONT,'fontweight',P.YLABELFONTWEIGHT,...
    'fontsize',P.YLABELFONTSIZE,'color',P.YLABELFONTCOLOR,...
    'interpreter',P.YLABELINTERP);

end
