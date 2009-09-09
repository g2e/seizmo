function [fh,sfh]=plot1(data,varargin)
%PLOT1    Plot SEIZMO data records in individual subplots
%
%    Usage:  [fh,sfh]=plot1(data,'plot_option',plot_option_value,...)
%
%    Description: PLOT1(DATA) plots timeseries and xy SEIZMO records in a
%     figure as individual subplots.  Other filetypes are ignored (a space 
%     is left for their subplot in the figure).  By default, markers 'o',
%     'a','f','t(n)' are drawn if defined and are given labels using 'ko',
%     'ka','kf','kt(n)' if defined otherwise the field name (o,a,f,t(n)) is
%     used.  
%
%     [FH,SFH]=PLOT1(DATA) returns the figure handle in FH and the subplot
%     handles in SFH.
%
%     PLOT1(DATA,PLOT_OPTION,OPTION_VALUE) modifies the plot parameter
%     PLOT_OPTION to OPTION_VALUE.  A list of available options and their
%     applicable range can be found with the function PCONF.
%
%    Examples:
%     Plot records with 3 columns of subplots:
%      plot1(data,'ncols',3)
%
%     Plot records with the 'jet' colormap:
%      plot1(data,'colormap','jet')
%
%     Create your own colormap (one color per row):
%      plot1(data,'colormap',['k'; 'r'; 'b'])
%       or
%      plot1(data,'colormap',[0 0 0; 1 0 0; 0 0 1])
%
%     Black on white plot without markers:
%      plot1(data,'bgcolor','w','fgcolor','k','colormap','k','markers',false)
%
%    See also:  plot2, plot0, recordsection

%     Version History:
%        ????????????? - Initial Version
%        June 12, 2008 - doc update
%        Apr. 23, 2009 - fix seizmocheck for octave, move usage up
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 23, 2009 at 20:35 GMT

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

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
        warning('seizmo:p1:badInput','Unknown Option: %s',varargin{i}); 
    end
end

% clean up unset parameters
P=plotconfigfix(P);

% select/open plot
if(isempty(P.FIGHANDLE) || P.FIGHANDLE<1)
    fh=figure;
else
    fh=figure(P.FIGHANDLE);
end

% SOME STYLING OF THE PLOT
set(gcf,'name',['P1 -- ' P.NAME],...
        'numbertitle',P.NUMBERTITLE,...
        'menubar',P.MENUBAR,...
        'toolbar',P.TOOLBAR,...
        'renderer',P.RENDERER,...
        'renderermode',P.RENDERERMODE,...
        'doublebuffer',P.DOUBLEBUFFER,...
        'color',P.FIGBGCOLOR,...
        'pointer',P.POINTER,...
        'units',P.UNITS,...
        'position',P.POSITION);

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
[t,kt,o,ko,a,ka,f,kf,b,e,delta,npts,gcarc]=...
    getheader(data,'t','kt','o','ko','a','ka','f','kf',...
    'b','e','delta','npts','gcarc');

% header structures (for determining if undefined)
[h,idx]=versioninfo(data);

% default columns/rows
if(isempty(P.NCOLS))
    P.NCOLS=fix(sqrt(nrecs));
    nrows=ceil(nrecs/P.NCOLS);
else
    nrows=ceil(nrecs/P.NCOLS);
end

% force x/y limits into individual subplot control
if(~isempty(P.XLIMITS))
    P.XLIMITS=repmat(P.XLIMITS,ceil(nrecs/size(P.XLIMITS,1)),1);
end
if(~isempty(P.YLIMITS))
    P.YLIMITS=repmat(P.XLIMITS,ceil(nrecs/size(P.XLIMITS,1)),1);
end

% loop through each record
sfh=zeros(nrecs,1);
for i=1:nrecs
    % header version
    v=idx(i);
    
    % check if timeseries or xy (if not skip)
    if(~any(strcmp(iftype(i),{'Time Series File' ...
            'General X vs Y file'})))
        continue; 
    end
    
    % get record timing
    if(strcmp(leven(i),'true'))
        time=b(i)+(0:npts(i)-1).'*delta(i);
    else
        time=data(i).ind; 
    end
    
    % focus to new subplot and draw record
    sfh(i)=subplot(nrows,P.NCOLS,i);
    plot(time,data(i).dep,'color',colors(i,:),'linewidth',P.RECWIDTH);
    set(gca,'fontname',P.AXISFONT,'fontweight',P.AXISFONTWEIGHT,...
        'fontsize',P.AXISFONTSIZE,'box',P.BOX,...
        'xcolor',P.XAXISCOLOR,'ycolor',P.YAXISCOLOR,...
        'xaxislocation',P.XAXISLOCATION,'yaxislocation',P.YAXISLOCATION,...
        'tickdir',P.TICKDIR,'ticklength',P.TICKLEN,...
        'linewidth',P.AXISLINEWIDTH,'color',P.PLOTBGCOLOR,...
        'xminortick',P.XMINORTICK,'yminortick',P.YMINORTICK);
    grid(P.GRID);
    
    % zooming
    axis(P.AXIS{:});
    if(~isempty(P.XLIMITS)); axis auto; xlim(P.XLIMITS(i,:)); end
    if(~isempty(P.YLIMITS)); ylim(P.YLIMITS(i,:)); end
    
    % scaling parameters for header flags
    fylimits=get(gca,'ylim');
    fxlimits=get(gca,'xlim');
    yrange=fylimits(2)-fylimits(1);
    xrange=fxlimits(2)-fxlimits(1);
    ypad=P.MARKYPAD*yrange;  % these are simplistic - really you will have to
    xpad=P.MARKXPAD*xrange;  % adjust these based on your font and plot size
    
    % axis labels
    if(isempty(P.TITLE))
        if(isfield(data,'name') && ~isempty(data(i).name))
            try
                p1title=[texlabel(data(i).name,'literal') ...
                    '  -  ' num2str(gcarc(i)) '\circ'];
            catch
    	        p1title=[data(i).name ...
        	    '  -  ' num2str(gcarc(i)) '\circ'];
            end
        else p1title=['RECORD ' num2str(i) ...
                '  -  ' num2str(gcarc(i)) '\circ'];
        end
    else
        p1title=P.TITLE;
    end
    if(isempty(P.XLABEL)); P.XLABEL='Time (sec)'; end
    if(isempty(P.YLABEL)); P.YLABEL='Amplitude'; end
    title(p1title,'fontname',P.TITLEFONT,'fontweight',P.TITLEFONTWEIGHT,...
        'fontsize',P.TITLEFONTSIZE,'color',P.TITLEFONTCOLOR,...
        'interpreter',P.TITLEINTERP);
    xlabel(P.XLABEL,'fontname',P.XLABELFONT,'fontweight',P.XLABELFONTWEIGHT,...
        'fontsize',P.XLABELFONTSIZE,'color',P.XLABELFONTCOLOR,...
        'interpreter',P.XLABELINTERP);
    ylabel(P.YLABEL,'fontname',P.YLABELFONT,'fontweight',P.YLABELFONTWEIGHT,...
        'fontsize',P.YLABELFONTSIZE,'color',P.YLABELFONTCOLOR,...
        'interpreter',P.YLABELINTERP);
    
    % skip markers
    if(~P.MARKERS); continue; end
    
    % plot origin flag
    hold on
    if(o(i)~=h(v).undef.ntype)
        plot([o(i) o(i)].',[fylimits(1)+ypad fylimits(2)-2*ypad].',...
            'color',P.OCOLOR,'linewidth',P.MARKERWIDTH);
        if(P.MARKERLABELS)
            if(~strcmp(ko{i},h(v).undef.stype))
                text(o(i)+xpad,fylimits(2)-2*ypad,ko{i},'color',P.OCOLOR, ...
                    'fontname',P.MARKERFONT,'fontweight',P.MARKERFONTWEIGHT,...
                    'fontsize',P.MARKERFONTSIZE,'color',P.OMARKERFONTCOLOR,...
                    'verticalalignment','top','clipping','on');
            else
                text(o(i)+xpad,fylimits(2)-2*ypad,'o','color',P.OCOLOR, ...
                    'fontname',P.MARKERFONT,'fontweight',P.MARKERFONTWEIGHT,...
                    'fontsize',P.MARKERFONTSIZE,'color',P.OMARKERFONTCOLOR,...
                    'verticalalignment','top','clipping','on');
            end
        end
    end
    
    % plot arrival flag
    if(a(i)~=h(v).undef.ntype)
        plot([a(i) a(i)].',[fylimits(1)+ypad fylimits(2)-2*ypad].',...
            'color',P.ACOLOR,'linewidth',P.MARKERWIDTH);
        if(P.MARKERLABELS)
            if(~strcmp(ka{i},h(v).undef.stype))
                text(a(i)+xpad,fylimits(2)-2*ypad,ka{i},'color',P.ACOLOR, ...
                    'fontname',P.MARKERFONT,'fontweight',P.MARKERFONTWEIGHT,...
                    'fontsize',P.MARKERFONTSIZE,'color',P.AMARKERFONTCOLOR,...
                    'verticalalignment','top','clipping','on');
            else
                text(a(i)+xpad,fylimits(2)-2*ypad,'a','color',P.ACOLOR, ...
                    'fontname',P.MARKERFONT,'fontweight',P.MARKERFONTWEIGHT,...
                    'fontsize',P.MARKERFONTSIZE,'color',P.AMARKERFONTCOLOR,...
                    'verticalalignment','top','clipping','on');
            end
        end
    end
    
    % plot finish flag
    if(f(i)~=h(v).undef.ntype)
        plot([f(i) f(i)].',[fylimits(1)+ypad fylimits(2)-2*ypad].',...
            'color',P.FCOLOR,'linewidth',P.MARKERWIDTH);
        if(P.MARKERLABELS)
            if(~strcmp(kf{i},h(v).undef.stype))
                text(f(i)+xpad,fylimits(2)-2*ypad,kf{i},'color',P.FCOLOR, ...
                    'fontname',P.MARKERFONT,'fontweight',P.MARKERFONTWEIGHT,...
                    'fontsize',P.MARKERFONTSIZE,'color',P.FMARKERFONTCOLOR,...
                    'verticalalignment','top','clipping','on');
            else
                text(f(i)+xpad,fylimits(2)-2*ypad,'f','color',P.FCOLOR, ...
                    'fontname',P.MARKERFONT,'fontweight',P.MARKERFONTWEIGHT,...
                    'fontsize',P.MARKERFONTSIZE,'color',P.FMARKERFONTCOLOR,...
                    'verticalalignment','top','clipping','on');
            end
        end
    end
    
    % plot picks
    for j=0:9
        if(t(i,j+1)~=h(v).undef.ntype)
            plot([t(i,j+1) t(i,j+1)].',...
                [fylimits(1)+ypad fylimits(2)-(1+mod(j,5))*ypad].',...
                'color',P.TCOLOR,'linewidth',P.MARKERWIDTH)
            if(P.MARKERLABELS)
                if(~strcmp(kt{i,j+1},h(v).undef.stype))
                    text(t(i,j+1)+xpad,...
                        fylimits(2)-(1+mod(j,5))*ypad,kt{i,j+1},'color',P.TCOLOR,...
                        'fontname',P.MARKERFONT,'fontweight',P.MARKERFONTWEIGHT,...
                        'fontsize',P.MARKERFONTSIZE,'color',P.TMARKERFONTCOLOR,...
                        'verticalalignment','top','clipping','on');
                else
                    text(t(i,j+1)+xpad,...
                        fylimits(2)-(1+mod(j,5))*ypad,['t' num2str(j)],'color',P.TCOLOR,...
                        'fontname',P.MARKERFONT,'fontweight',P.MARKERFONTWEIGHT,...
                        'fontsize',P.MARKERFONTSIZE,'color',P.TMARKERFONTCOLOR,...
                        'verticalalignment','top','clipping','on');
                end
            end
        end
    end

    % rezooming
    axis(P.AXIS{:});
    xlim([min(b) max(e)]);
    if(~isempty(P.XLIMITS)); axis auto; xlim(P.XLIMITS(i,:)); end
    if(~isempty(P.YLIMITS)); ylim(P.YLIMITS(i,:)); end
    hold off
end

end
