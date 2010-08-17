function [varargout]=plot1(data,varargin)
%PLOT1    Plot SEIZMO data records in individual subplots
%
%    Usage:    plot1(data)
%              plot1(...,'option',value,...)
%              ax=plot1(...)
%
%    Description:
%     PLOT1(DATA) draws all non-xyz records in SEIZMO struct DATA in
%     separate subplots in a new figure.  Each record is plotted as a
%     distinct color from the HSV colormap.  Spectral records are converted
%     to the time domain prior to plotting.
%
%     PLOT1(...,'OPTION',VALUE,...) sets certain plotting options to do
%     simple manipulation of the plots.  Available options are:
%      FGCOLOR    -- foreground color (axes, text, labels)
%      BGCOLOR    -- background color (does not set figure color)
%      AXIS       -- axes to plot in
%      COLORMAP   -- colormap for coloring data
%      XLABEL     -- x axis label
%      YLABEL     -- y axis label
%      TITLE      -- title
%      XLIM       -- x axis limits (tight by default)
%      YLIM       -- y axis limits (tight by default)
%      LINEWIDTH  -- line width of records (default is 1)
%      LINESTYLE  -- line style of records (can be char/cellstr array)
%      NUMCOLS    -- number of subplot columns
%      UTC        -- plot in absolute time if TRUE (UTC, no leap support)
%      DATEFORMAT -- date format used if ABSOLUTE (default is auto)
%      XDIR       -- 'normal' or 'reverse'
%      YDIR       -- 'normal' or 'reverse'
%      FONTSIZE   -- size of fonts in the axes
%      ALIGN      -- ignore label and tick overlaps when aligning subplots
%      XSCALE     -- 'linear' or 'log'
%      YSCALE     -- 'linear' or 'log'
%
%     AX=PLOT1(...) returns the handles for all the axes drawn in.  This is
%     useful for more detailed plot manipulation.
%
%    Notes:
%
%    Examples:
%     % plot records with 3 columns of subplots:
%     plot1(data,'ncols',3)
%
%     % plot records with the 'jet' colormap:
%     plot1(data,'colormap','jet')
%
%     % create your own colormap (one color per row):
%     plot1(data,'colormap',['k'; 'r'; 'b'])
%       or
%     plot1(data,'colormap',[0 0 0; 1 0 0; 0 0 1])
%
%    See also: PLOT0, PLOT2, RECORDSECTION

%     Version History:
%        Aug. 14, 2010 - rewrite
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 14, 2010 at 23:00 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% check struct
error(seizmocheck(data,'dep'));
nrecs=numel(data);

% default/parse options
opt=parse_seizmo_plot_options(varargin{:});

% line coloring
if(ischar(opt.CMAP) || iscellstr(opt.CMAP))
    % list of color names or a colormap function
    try
        % attempt color name to rgb conversion first
        opt.CMAP=name2rgb(opt.CMAP);
        opt.CMAP=repmat(opt.CMAP,ceil(nrecs/size(opt.CMAP,1)),1);
    catch
        % guess its a colormap function then
        opt.CMAP=str2func(opt.CMAP);
        opt.CMAP=opt.CMAP(nrecs);
    end
else
    % numeric colormap array
    opt.CMAP=repmat(opt.CMAP,ceil(nrecs/size(opt.CMAP,1)),1);
end

% line style
opt.LINESTYLE=cellstr(opt.LINESTYLE);
opt.LINESTYLE=opt.LINESTYLE(:);
opt.LINESTYLE=repmat(opt.LINESTYLE,ceil(nrecs/size(opt.LINESTYLE,1)),1);

% check filetype (only timeseries or xy)
iftype=getenumid(data,'iftype');
time=strcmpi(iftype,'itime') | strcmpi(iftype,'ixy');
spec=strcmpi(iftype,'irlim') | strcmpi(iftype,'iamph');
goodfiles=find(time | spec)';

% convert spectral to timeseries
if(sum(spec)); data(spec)=idft(data(spec)); end

% header info
leven=getlgc(data,'leven');
idep=shortidep(getenumdesc(data,'idep'));
[b,npts,delta,z6,kname]=getheader(data,...
    'b','npts','delta','z6','kname');
z6=datenum(cell2mat(z6));

% subplot style
if(isempty(opt.AXIS) || numel(opt.AXIS)~=nrecs || ~isreal(opt.AXIS) ...
        || any(~ishandle(opt.AXIS)) ...
        || any(~strcmp('axes',get(opt.AXIS,'type'))))
    % new figure
    fh=figure('color',opt.BGCOLOR);
    if(isempty(opt.NUMCOLS))
        opt.NUMCOLS=fix(sqrt(nrecs));
    end
    nrows=ceil(nrecs/opt.NUMCOLS);
    opt.AXIS=makesubplots(nrows,opt.NUMCOLS,1:nrecs,...
        opt.ALIGN{:},'parent',fh);
end

% loop through every record
for i=goodfiles
    % adjust current axis
    cla(opt.AXIS(i),'reset');
    set(opt.AXIS(i),'ydir',opt.YDIR,'xdir',opt.XDIR,...
        'xscale',opt.XSCALE,'yscale',opt.YSCALE,...
        'fontsize',opt.FONTSIZE,'color',opt.BGCOLOR,...
        'xcolor',opt.FGCOLOR,'ycolor',opt.FGCOLOR);
    
    % plot record
    hold(opt.AXIS(i),'on');
    switch leven{i}
        case 'false'
            if(opt.ABSOLUTE)
                plot(opt.AXIS(i),z6(i)+data(i).ind/86400,data(i).dep,...
                    'color',opt.CMAP(i,:),...
                    'linestyle',opt.LINESTYLE{i},...
                    'linewidth',opt.LINEWIDTH);
            else
                plot(opt.AXIS(i),data(i).ind,data(i).dep,...
                    'color',opt.CMAP(i,:),...
                    'linestyle',opt.LINESTYLE{i},...
                    'linewidth',opt.LINEWIDTH);
            end
        otherwise
            if(opt.ABSOLUTE)
                plot(opt.AXIS(i),...
                    z6(i)+(b(i)+(0:npts(i)-1)*delta(i)).'/86400,...
                    data(i).dep,...
                    'color',opt.CMAP(i,:),...
                    'linestyle',opt.LINESTYLE{i},...
                    'linewidth',opt.LINEWIDTH);
            else
                plot(opt.AXIS(i),(b(i)+(0:npts(i)-1)*delta(i)).',...
                    data(i).dep,...
                    'color',opt.CMAP(i,:),...
                    'linestyle',opt.LINESTYLE{i},...
                    'linewidth',opt.LINEWIDTH);
            end
    end
    hold(opt.AXIS(i),'off');
    
    % extras
    box(opt.AXIS(i),'on');
    grid(opt.AXIS(i),'on');
    axis(opt.AXIS(i),'tight');
    
    % axis zooming
    if(~isempty(opt.XLIM))
        %if(isempty(opt.YLIM)); axis autoy; end
        xlim(opt.AXIS(i),opt.XLIM);
    end
    if(~isempty(opt.YLIM))
        %if(isempty(opt.XLIM)); axis autox; end
        ylim(opt.AXIS(i),opt.YLIM);
    end
    
    % datetick
    if(opt.ABSOLUTE)
        if(isempty(opt.DATEFORMAT))
            if(isempty(opt.XLIM))
                datetick(opt.AXIS(i),'x');
            else
                datetick(opt.AXIS(i),'x','keeplimits');
            end
        else
            if(isempty(opt.XLIM))
                datetick(opt.AXIS(i),'x',opt.DATEFORMAT);
            else
                datetick(opt.AXIS(i),'x',opt.DATEFORMAT,'keeplimits');
            end
        end
    end
    
    % label
    if(isempty(opt.TITLE))
        if(~isempty(data(i).name))
            p1title=texlabel(data(i).name,'literal');
        else
            p1title=['RECORD ' num2str(i)];
        end
    elseif(isnumeric(opt.TITLE))
        switch opt.TITLE
            case 1 % kstnm
                p1title=kname(i,2);
            case 2 % stcmp
                p1title=strcat(kname(i,2),'.',kname(i,4));
            case 3 % kname
                p1title=texlabel(strcat(kname(i,1),'.',kname(i,2),...
                '.',kname(i,3),'.',kname(i,4)),'literal');
            otherwise
                p1title=['RECORD ' num2str(i)];
        end
    else
        p1title=opt.TITLE;
    end
    if(isempty(opt.XLABEL))
        if(opt.ABSOLUTE)
            xlimits=get(opt.AXIS(i),'xlim');
            p1xlabel=joinwords(cellstr(datestr(unique(fix(xlimits)))),...
                '   to   ');
        else
            p1xlabel='Time (sec)';
        end
    else
        p1xlabel=opt.XLABEL;
    end
    if(isempty(opt.YLABEL))
        p1ylabel=idep(i);
    else
        p1ylabel=opt.YLABEL;
    end
    title(opt.AXIS(i),p1title,...
        'color',opt.FGCOLOR,'fontsize',opt.FONTSIZE);
    xlabel(opt.AXIS(i),p1xlabel,...
        'color',opt.FGCOLOR,'fontsize',opt.FONTSIZE);
    ylabel(opt.AXIS(i),p1ylabel,...
        'color',opt.FGCOLOR,'fontsize',opt.FONTSIZE);
end

% output axes if wanted
if(nargout); varargout{1}=opt.AXIS; end

end


function [idep]=shortidep(idep)
%SHORTIDEP    Converts idep long description to just units
long={'Unknown' 'Displacement (nm)' 'Velocity (nm/sec)' ...
      'Acceleration (nm/sec^2)' 'Velocity (volts)' 'Absement (nm*sec)' ...
      'Absity (nm*sec^2)' 'Abseleration (nm*sec^3)' 'Abserk (nm*sec^4)' ...
      'Absnap (nm*sec^5)' 'Absackle (nm*sec^6)' 'Abspop (nm*sec^7)' ...
      'Jerk (nm/sec^3)' 'Snap (nm/sec^4)' 'Crackle (nm/sec^5)' ...
      'Pop (nm/sec^6)' 'Counts'};
short={'unknown' 'nm' 'nm/sec' 'nm/sec^2' 'volts' 'nm*sec' 'nm*sec^2' ...
       'nm*sec^3' 'nm*sec^4' 'nm*sec^5' 'nm*sec^6' 'nm*sec^7' ...
       'nm/sec^3' 'nm/sec^4' 'nm/sec^5' 'nm/sec^6' 'counts'};
[tf,i]=ismember(idep,long);
idep(~tf)={'unknown'};
idep(tf)=short(i(tf));
end

