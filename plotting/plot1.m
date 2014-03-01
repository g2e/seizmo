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
%      FONTWEIGHT -- 'light', 'normal', 'demi' or 'bold'
%      ALIGN      -- ignore label and tick overlaps when aligning subplots
%      XSCALE     -- 'linear' or 'log'
%      YSCALE     -- 'linear' or 'log'
%      MARKERS    -- true/false where true draws the markers
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
%        Sep. 14, 2010 - added marker support
%        Feb. 10, 2011 - clear todo list (was about markers)
%        Apr. 16, 2011 - allow empty title/xlabel/ylabel, allow datetick
%                        for non-absolute
%        Apr. 19, 2011 - userdata for each record contains record metadata
%        Nov.  8, 2011 - move shortidep out as an independent function
%        Nov. 11, 2011 - per record linewidth, fontweight support
%        Feb.  6, 2012 - better getheader usage, set displayname to kname
%                        string, fix legend coloring
%        Aug. 30, 2012 - allow [] input for labeling, short idep by default
%        Mar.  1, 2014 - maketime fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 23:00 GMT

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

% line style/width
opt.LINESTYLE=cellstr(opt.LINESTYLE);
opt.LINESTYLE=opt.LINESTYLE(:);
opt.LINESTYLE=repmat(opt.LINESTYLE,ceil(nrecs/size(opt.LINESTYLE,1)),1);
opt.LINEWIDTH=opt.LINEWIDTH(:);
opt.LINEWIDTH=repmat(opt.LINEWIDTH,ceil(nrecs/size(opt.LINEWIDTH,1)),1);

% check filetype (only timeseries or xy)
iftype=getheader(data,'iftype id');
time=strcmpi(iftype,'itime') | strcmpi(iftype,'ixy');
spec=strcmpi(iftype,'irlim') | strcmpi(iftype,'iamph');
goodfiles=find(time | spec)';

% convert spectral to timeseries
if(sum(spec)); data(spec)=idft(data(spec)); end

% header info
[b,npts,delta,z6,kname,leven,idep]=getheader(data,...
    'b','npts','delta','z6','kname','leven lgc','idep desc');
z6=datenum(cell2mat(z6));

% convert idep to short form
idep=shortidep(idep);

% names for legend
displayname=strcat(kname(:,1),'.',kname(:,2),...
    '.',kname(:,3),'.',kname(:,4));

% get markers info
[marknames,marktimes]=getmarkers(data);

% convert markers to absolute time if used
if(opt.ABSOLUTE)
    marktimes=marktimes/86400+z6(:,ones(1,size(marktimes,2)));
end

 % subplot style
if(isempty(opt.AXIS) || numel(opt.AXIS)~=nrecs || ~isreal(opt.AXIS) ...
        || any(~ishandle(opt.AXIS)) ...
        || any(~strcmp('axes',get(opt.AXIS,'type'))))
    % new figure
    fh=figure('color',opt.BGCOLOR,'defaulttextcolor',opt.FGCOLOR,...
        'defaultaxesxcolor',opt.FGCOLOR); % defaults for legend
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
        'fontsize',opt.FONTSIZE,'fontweight',opt.FONTWEIGHT,...
        'color',opt.BGCOLOR,'xcolor',opt.FGCOLOR,'ycolor',opt.FGCOLOR);
    
    % plot record
    hold(opt.AXIS(i),'on');
    switch leven{i}
        case 'false'
            if(opt.ABSOLUTE)
                rh=plot(opt.AXIS(i),z6(i)+data(i).ind/86400,data(i).dep,...
                    'color',opt.CMAP(i,:),...
                    'linestyle',opt.LINESTYLE{i},...
                    'linewidth',opt.LINEWIDTH(i));
            else
                rh=plot(opt.AXIS(i),data(i).ind,data(i).dep,...
                    'color',opt.CMAP(i,:),...
                    'linestyle',opt.LINESTYLE{i},...
                    'linewidth',opt.LINEWIDTH(i));
            end
        otherwise
            if(opt.ABSOLUTE)
                rh=plot(opt.AXIS(i),...
                    z6(i)+b(i)/86400+...
                    (0:delta(i)/86400:delta(i)/86400*(npts(i)-1)).',...
                    data(i).dep,...
                    'color',opt.CMAP(i,:),...
                    'linestyle',opt.LINESTYLE{i},...
                    'linewidth',opt.LINEWIDTH(i));
            else
                rh=plot(opt.AXIS(i),b(i)+...
                    (0:delta(i):delta(i)*(npts(i)-1)).',...
                    data(i).dep,...
                    'color',opt.CMAP(i,:),...
                    'linestyle',opt.LINESTYLE{i},...
                    'linewidth',opt.LINEWIDTH(i));
            end
    end
    hold(opt.AXIS(i),'off');
    
    % set userdata to everything but data (cleared)
    data(i).dep=[];
    data(i).ind=[];
    nrh=numel(rh);
    for ridx=1:nrh
        if(nrh>1)
            ncmpstr=[' (' num2str(ridx) ')'];
        else
            ncmpstr=[];
        end
        data(i).index=[i ridx];
        set(rh(ridx),'userdata',data(i),...
            'displayname',[displayname{i} ncmpstr]);
    end
    
    % tag records
    set(rh,'tag','record');
    
    % extras
    box(opt.AXIS(i),'on');
    grid(opt.AXIS(i),'on');
    axis(opt.AXIS(i),'tight');
    
    % add markers to axis userdata
    userdata.markers.names=marknames(i,:);
    userdata.markers.times=marktimes(i,:);
    userdata.function='plot1';
    set(opt.AXIS(i),'userdata',userdata);
    
    % draw markers
    if(opt.MARKERS)
        drawmarkers(opt.AXIS(i),varargin{:});
    end
    
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
    else
        if(~isempty(opt.DATEFORMAT))
            if(isempty(opt.XLIM))
                datetick(opt.AXIS(i),'x',opt.DATEFORMAT);
            else
                datetick(opt.AXIS(i),'x',opt.DATEFORMAT,'keeplimits');
            end
        end
    end
    
    % label
    if(~isempty(opt.TITLE) && isnumeric(opt.TITLE))
        switch opt.TITLE
            case 1 % filename
                if(~isempty(data(i).name))
                    if(isoctave)
                        p1title=data(i).name;
                    else
                        p1title=texlabel(data(i).name,'literal');
                    end
                else
                    p1title=['RECORD ' num2str(i)];
                end
            case 2 % kstnm
                p1title=kname(i,2);
            case 3 % kcmpnm
                p1title=kname(i,4);
            case 4 % shortcmp
                p1title=kname{i,4}(3);
            case 5 % stashort
                p1title=strcat(kname(i,2),'.',kname{i,4}(3));
            case 6 % stcmp
                p1title=strcat(kname(i,2),'.',kname(i,4));
            case 7 % kname
                p1title=texlabel(strcat(kname(i,1),'.',kname(i,2),...
                '.',kname(i,3),'.',kname(i,4)),'literal');
            otherwise
                p1title=['RECORD ' num2str(i)];
        end
    else
        p1title=opt.TITLE;
    end
    if(~isempty(opt.XLABEL) && isnumeric(opt.XLABEL) && opt.XLABEL==1)
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
    if(~isempty(opt.YLABEL) && isnumeric(opt.YLABEL) && opt.YLABEL==1)
        p1ylabel=idep(i);
    else
        p1ylabel=opt.YLABEL;
    end
    title(opt.AXIS(i),p1title,'color',opt.FGCOLOR,...
        'fontsize',opt.FONTSIZE,'fontweight',opt.FONTWEIGHT);
    xlabel(opt.AXIS(i),p1xlabel,'color',opt.FGCOLOR,...
        'fontsize',opt.FONTSIZE,'fontweight',opt.FONTWEIGHT);
    ylabel(opt.AXIS(i),p1ylabel,'color',opt.FGCOLOR,...
        'fontsize',opt.FONTSIZE,'fontweight',opt.FONTWEIGHT);
end

% output axes if wanted
if(nargout); varargout{1}=opt.AXIS; end

end
