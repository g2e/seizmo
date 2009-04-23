function [fh,lh]=plot0(data,varargin)
%PLOT0   Plot SEIZMO data records in an evenly spaced record section
%
%    Usage:  [fh,lh]=plot0(data,'plot_option',plot_option_value,...)
%
%    Description: Plots timeseries and xy SEIZMO records spaced out evenly
%     in a single plot.  Other record types are ignored.  Optional
%     inputs should correspond to fields returned by function pconf.
%     Outputs are the figure and legend handles.
%
%    Examples:
%     To add record names to the yaxis:
%        plot0(data,'namesonyaxis',true)
%
%    See also:  plot1, plot2, recordsection

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
        warning('seizmo:p3:badInput','Unknown Option: %s',varargin{i}); 
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
set(gcf,'name',['P3 -- ' P.NAME],...
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
[t,kt,o,ko,a,ka,f,kf,b,e,npts,delta,depmin,depmax]=...
    getheader(data,'t','kt','o','ko','a','ka','f','kf',...
    'b','e','npts','delta','depmin','depmax');

% header structures (for determining if undefined)
[h,idx]=versioninfo(data);

% yaxis scaling for amplitudes
if(P.NORM2YRANGE)
    scale=nrecs*P.NORMMAX;
else
    scale=P.NORMMAX;
end

% check normalization style
if(~any(strcmpi(P.NORMSTYLE,{'single' 'group'})))
    warning('seizmo:p3:badInput','bad normalization style')
    P.NORMSTYLE='single';
end

% check filetype (only timeseries or xy)
goodfiles=(strcmp(iftype,'Time Series File') | ...
    strcmp(iftype,'General X vs Y file'));
indices=find(goodfiles).';

% find amplitude scaling
ampmax=max(abs([depmin.'; depmax.';]));
if(strcmpi(P.NORMSTYLE,'group')); ampmax(:)=max(ampmax(indices)); end
ampmax(ampmax==0)=1; % avoid NaNs

% find time scaling
xpad=0.005*abs(max(e)-min(b));

% loop through each file
hold on
for i=indices
    % header version
    v=idx(i);
    
    % get record timing
    if(strcmp(leven(i),'true')); time=(b(i)+(0:npts(i)-1)*delta(i)).';
    else time=data(i).ind; end
    
    % plot series
    plot(time,i+data(i).dep/ampmax(i)*scale,...
        'color',colors(i,:),'linewidth',P.RECWIDTH);
    
    % skip markers
    if(~P.MARKERS); continue; end
    
    % plot origin flag
    hold on
    if(o(i)~=h(v).undef.ntype)
        plot([o(i) o(i)].',i+[-0.4 0.4].',...
            'color',P.OCOLOR,'linewidth',P.MARKERWIDTH);
        if(P.MARKERLABELS)
            if(~strcmp(ko{i},h(v).undef.stype))
                text(o(i)+xpad,i+0.3,ko{i},'color',P.OCOLOR, ...
                    'fontname',P.MARKERFONT,'fontweight',P.MARKERFONTWEIGHT,...
                    'fontsize',P.MARKERFONTSIZE,'color',P.OMARKERFONTCOLOR,...
                    'verticalalignment','top','clipping','on');
            else
                text(o(i)+xpad,i+0.3,'o','color',P.OCOLOR, ...
                    'fontname',P.MARKERFONT,'fontweight',P.MARKERFONTWEIGHT,...
                    'fontsize',P.MARKERFONTSIZE,'color',P.OMARKERFONTCOLOR,...
                    'verticalalignment','top','clipping','on');
            end
        end
    end
    
    % plot arrival flag
    if(a(i)~=h(v).undef.ntype)
        plot([a(i) a(i)].',i+[-0.4 0.4].',...
            'color',P.ACOLOR,'linewidth',P.MARKERWIDTH);
        if(P.MARKERLABELS)
            if(~strcmp(ka{i},h(v).undef.stype))
                text(a(i)+xpad,i+0.3,ka{i},'color',P.ACOLOR, ...
                    'fontname',P.MARKERFONT,'fontweight',P.MARKERFONTWEIGHT,...
                    'fontsize',P.MARKERFONTSIZE,'color',P.AMARKERFONTCOLOR,...
                    'verticalalignment','top','clipping','on');
            else
                text(a(i)+xpad,i+0.3,'a','color',P.ACOLOR, ...
                    'fontname',P.MARKERFONT,'fontweight',P.MARKERFONTWEIGHT,...
                    'fontsize',P.MARKERFONTSIZE,'color',P.AMARKERFONTCOLOR,...
                    'verticalalignment','top','clipping','on');
            end
        end
    end
    
    % plot finish flag
    if(f(i)~=h(v).undef.ntype)
        plot([f(i) f(i)].',i+[-0.4 0.4].',...
            'color',P.FCOLOR,'linewidth',P.MARKERWIDTH);
        if(P.MARKERLABELS)
            if(~strcmp(kf{i},h(v).undef.stype))
                text(f(i)+xpad,i+0.3,kf{i},'color',P.FCOLOR, ...
                    'fontname',P.MARKERFONT,'fontweight',P.MARKERFONTWEIGHT,...
                    'fontsize',P.MARKERFONTSIZE,'color',P.FMARKERFONTCOLOR,...
                    'verticalalignment','top','clipping','on');
            else
                text(f(i)+xpad,i+0.3,'f','color',P.FCOLOR, ...
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
                i+[-0.4 0.4].',...
                'color',P.TCOLOR,'linewidth',P.MARKERWIDTH)
            if(P.MARKERLABELS)
                if(~strcmp(kt{i,j+1},h(v).undef.stype))
                    text(t(i,j+1)+xpad,...
                        i+0.3,kt{i,j+1},'color',P.TCOLOR,...
                        'fontname',P.MARKERFONT,'fontweight',P.MARKERFONTWEIGHT,...
                        'fontsize',P.MARKERFONTSIZE,'color',P.TMARKERFONTCOLOR,...
                        'verticalalignment','top','clipping','on');
                else
                    text(t(i,j+1)+xpad,...
                        i+0.3,['t' num2str(j)],'color',P.TCOLOR,...
                        'fontname',P.MARKERFONT,'fontweight',P.MARKERFONTWEIGHT,...
                        'fontsize',P.MARKERFONTSIZE,'color',P.TMARKERFONTCOLOR,...
                        'verticalalignment','top','clipping','on');
                end
            end
        end
    end
end
hold off

% zooming
axis(P.AXIS{:});
if(~isempty(P.XLIMITS)); axis auto; xlim(P.XLIMITS); end
if(~isempty(P.YLIMITS)); ylim(P.YLIMITS); end

% label y-axis with record names
if(P.NAMESONYAXIS)
    if(isfield(data,'name'))
        set(gca,'ytick',1:nrecs);
        set(gca,'yticklabel',{data.name}.');
    end
end

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
if(isempty(P.YLABEL) && ~P.NAMESONYAXIS); P.YLABEL='Record Number'; end
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
