function [fh,lh]=p2(data,varargin)
%P2    Overlay plot of SEIZMO data records
%
%    Description: Plots timeseries and xy SEIZMO records over one another
%     in a single plot.  Other record types are ignored.  Optional
%     inputs should correspond to fields returned by function pconf.
%     Outputs are the figure and legend handles.
%
%    Usage:  [fh,lh]=p2(data,'plot_option',plot_option_value,...)
%
%    Examples:
%     To overlay the first 4 records
%      p2(data(1:4)) 
%
%     To overlay the 5th and 8th records from 0 to 300 seconds
%      p2(data([5 8]),'xlimits',[0 300]) 
%
%     To plot all traces with a legend, without limiting the x/y axis
%      p2(data,'legend',true)
%
%    See also:  p1, p3, recsec

% check data structure
error(seischk(data,'dep'))

% get plotting style defaults
P=pconf;

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
        warning('seizmo:p2:badInput','Unknown Option: %s',varargin{i}); 
    end
end

% clean up unset parameters
P=pconffix(P);

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
set(gcf,'name',['P2 -- ' P.NAME],...
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
leven=glgc(data,'leven');
error(lgcchk('leven',leven))
iftype=genumdesc(data,'iftype');
[b,npts,delta,depmin,depmax]=...
    gh(data,'b','npts','delta','depmin','depmax');

% check normalization style
if(~any(strcmpi(P.NORMSTYLE,{'single' 'group'})))
    warning('seizmo:p2:badInput','bad normalization style')
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
    plot(time,data(i).dep/(ampmax(i)^P.P2NORM)*(P.NORMMAX^P.P2NORM),...
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
if(isempty(P.YLABEL)); P.YLABEL='Amplitude'; end
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
