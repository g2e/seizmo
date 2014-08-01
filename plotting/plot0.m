function [varargout]=plot0(data,varargin)
%PLOT0    Evenly spaced plot of SEIZMO records
%
%    Usage:    plot0(data)
%              plot0(...,'option',value,...)
%              ax=plot0(...)
%
%    Description:
%     PLOT0(DATA) draws all non-xyz records in SEIZMO struct DATA in a new
%     figure in the same axes.  The records are normalized as a group (with
%     a maximum amplitude range of 1/3 the yaxis range) and spaced at unit
%     distance from one another.  Each record is drawn as a distinct color
%     from the HSV colormap.  Spectral records are converted to the time
%     domain prior to plotting.
%
%     PLOT0(...,'OPTION',VALUE,...) sets certain plotting options to do
%     simple manipulation of the plots.  Available options are:
%      FGCOLOR      -- foreground color (axes, text, labels)
%      BGCOLOR      -- background color (does not set figure color)
%      AXIS         -- axes to plot in
%      COLORMAP     -- colormap for coloring data
%      XLABEL       -- x axis label
%      YLABEL       -- y axis label
%      TITLE        -- title
%      XLIM         -- x axis limits (tight by default)
%      YLIM         -- y axis limits (tight by default)
%      LINEWIDTH    -- line width of records (default is 1)
%      LINESTYLE    -- line style of records (can be char/cellstr array)
%      NUMCOLS      -- number of subplot columns
%      UTC          -- plot in absolute time if TRUE (UTC, no leap support)
%      DATEFORMAT   -- date format used if ABSOLUTE (default is auto)
%      NORMSTYLE    -- normalize 'individually' or as a 'group'
%      NORMMAX      -- max value of normalized records
%      NORM2YAXIS   -- scale to yaxis range (NORMMAX is fraction of range)
%      NAMESONYAXIS -- true/false or 'kstnm' 'stcmp' 'kname'
%      XDIR         -- 'normal' or 'reverse'
%      YDIR         -- 'normal' or 'reverse'
%      FONTSIZE     -- size of fonts in the axes
%      FONTWEIGHT   -- 'light', 'normal', 'demi' or 'bold'
%      XSCALE       -- 'linear' or 'log'
%      YSCALE       -- 'linear' or 'log'
%      AMPSCALE     -- 'linear' or 'log'
%      MARKERS      -- true/false where true draws the markers
%
%     AX=PLOT0(...) returns the handle for the axis drawn in.  This is
%     useful for more detailed plot manipulation.
%
%    Notes:
%     - Using AMPSCALE='log' with negative data will find the absolute
%       value on the data.  If you try to recover the data with
%       EXTRACT_PLOT_DATA you will get all positive data.
%
%    Examples:
%     % add station+component names to the yaxis
%     plot0(data,'namesonyaxis','stcmp')
%
%    See also: PLOT1, PLOT2, RECORDSECTION

%     Version History:
%        Aug. 14, 2010 - rewrite
%        Sep. 14, 2010 - added marker support
%        Feb.  6, 2011 - fixed new figure/axes call
%        Apr. 16, 2011 - allow empty title/xlabel/ylabel, added several new
%                        nameonyaxis types, allow datetick for non-absolute
%        Apr. 19, 2011 - userdata for each record contains record metadata,
%                        drop the confusing upside-down default
%        Jan. 12, 2012 - minor improvement to normstyle handling
%        Jan. 25, 2012 - norm2yaxis bugfix
%        Feb.  6, 2012 - better getheader usage, set displayname to kname
%                        string, fix legend coloring, fix ncmp>1 brokeness
%        Mar.  6, 2012 - don't use depmin/depmax from header
%        Mar. 13, 2012 - make sure dep* is updated for extract_plot_data
%        Aug. 30, 2012 - allow [] input for labeling
%        Mar.  1, 2014 - maketime fix
%        July 25, 2014 - bugfix: handling of empty records
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 25, 2014 at 23:00 GMT

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

% check filetype (only timeseries or spectral)
iftype=getheader(data,'iftype id');
time=strcmpi(iftype,'itime') | strcmpi(iftype,'ixy');
spec=strcmpi(iftype,'irlim') | strcmpi(iftype,'iamph');
goodfiles=find(time | spec)';

% convert spectral to timeseries
if(sum(spec)); data(spec)=idft(data(spec)); end

% fix dep* for extraction later          
oldcheckheader_state=checkheader_state(true);
data=checkheader(data,'all','ignore','old_dep_stats','fix');
checkheader_state(oldcheckheader_state);

% header info
[b,npts,delta,z6,kname,leven,depmin,depmax]=getheader(data,...
    'b','npts','delta','z6','kname','leven lgc','depmin','depmax');
z6=datenum(cell2mat(z6));  % causes issues for record starts in leapsecond

% names for legend
displayname=strcat(kname(:,1),'.',kname(:,2),...
    '.',kname(:,3),'.',kname(:,4));

% get markers info
[marknames,marktimes]=getmarkers(data);

% convert markers to absolute time if used
if(opt.ABSOLUTE)
    marktimes=marktimes/86400+z6(:,ones(1,size(marktimes,2)));
end

% get normalizers based on data limits
switch lower(opt.AMPSCALE)
    case 'linear'
        depmin=abs(depmin);
        depmax=abs(depmax);
    case 'log'
        if(any(depmin<0))
            warning('seizmo:plot0:logOfNonPositive',...
                'Attempted logarithmic plotting of negative values!');
            for i=goodfiles
                if(depmin(i)<=0)
                    if(~isempty(data(i).dep))
                        depmin(i)=min(real(log10(...
                            data(i).dep(data(i).dep~=0))));
                        depmax(i)=max(real(log10(...
                            data(i).dep(data(i).dep~=0))));
                    end
                else
                    depmin(i)=log10(depmin(i));
                    depmax(i)=log10(depmax(i));
                end
            end
        else
            depmin=log10(depmin);
            depmax=log10(depmax);
        end
end

% find the vertical range of each record on the plot
if(opt.NORM2YAXIS)
    scale=nrecs*opt.NORMMAX/2;
else
    scale=opt.NORMMAX;
end
yrng=[(1:nrecs)' (1:nrecs)'];
yrng(:,1)=yrng(:,1)-scale;
yrng(:,2)=yrng(:,2)+scale;

% normalize
switch opt.NORMSTYLE
    case {1 'i' 'single' 'individually' 'individual' 'one' 'separately'}
        switch lower(opt.AMPSCALE)
            case 'linear'
                ampmax=max(depmin,depmax);
                ampmax(ampmax==0)=1; % avoid divide-by-zero below
                for i=goodfiles
                    data(i).dep=i+data(i).dep/ampmax(i)*scale;
                end
            case 'log'
                logrng=depmax-depmin;
                logrng(logrng==0)=1; % avoid divide-by-zero below
                for i=goodfiles
                    data(i).dep=i+(2*((real(log10(data(i).dep))...
                        -depmin(i))/logrng(i))-1)*scale;
                end
        end
    case {0 'g' 'a' 'group' 'together' 'all'}
        switch lower(opt.AMPSCALE)
            case 'linear'
                ampmax=max([depmin; depmax]);
                ampmax(ampmax==0)=1; % avoid divide-by-zero below
                for i=goodfiles
                    data(i).dep=i+data(i).dep/ampmax*scale;
                end
            case 'log'
                logmin=min(depmin);
                logmax=max(depmax);
                logrng=logmax-logmin;
                logrng(logrng==0)=1; % avoid divide-by-zero below
                for i=goodfiles
                    data(i).dep=i+(2*((real(log10(data(i).dep))...
                        -logmin)/logrng)-1)*scale;
                end
        end
    otherwise
        error('seizmo:plot0:badInput',...
            'Unknown NORMSTYLE value!');
end

% all in one plot
if(isempty(opt.AXIS) || ~isscalar(opt.AXIS) || ~isreal(opt.AXIS) ...
        || ~ishandle(opt.AXIS) || ~strcmp('axes',get(opt.AXIS,'type')))
    % new figure
    fh=figure('color',opt.BGCOLOR,'defaulttextcolor',opt.FGCOLOR,...
        'defaultaxesxcolor',opt.FGCOLOR); % defaults for legend
    opt.AXIS=axes('parent',fh);
else
    cla(opt.AXIS,'reset');
end

% adjust current axis
set(opt.AXIS,'ydir',opt.YDIR,'xdir',opt.XDIR,...
    'xscale',opt.XSCALE,'yscale',opt.YSCALE,...
    'fontsize',opt.FONTSIZE,'color',opt.BGCOLOR,...
    'xcolor',opt.FGCOLOR,'ycolor',opt.FGCOLOR);

% setup for markers
userdata.markers.records=nan(nrecs,1);

% for extraction
[data.empty]=deal(false);

% loop through every record
hold(opt.AXIS,'on');
for i=1:nrecs
    % handle empty and bad
    if(isempty(data(i).dep) || ~goodfiles(i))
        rh=plot(opt.AXIS,nan,i,...
            'color',opt.CMAP(i,:),...
            'linestyle',opt.LINESTYLE{i},...
            'linewidth',opt.LINEWIDTH);
    else
        switch leven{i}
            case 'false'
                if(opt.ABSOLUTE)
                    rh=plot(opt.AXIS,...
                        z6(i)+data(i).ind/86400,data(i).dep,...
                        'color',opt.CMAP(i,:),...
                        'linestyle',opt.LINESTYLE{i},...
                        'linewidth',opt.LINEWIDTH);
                else
                    rh=plot(opt.AXIS,data(i).ind,data(i).dep,...
                        'color',opt.CMAP(i,:),...
                        'linestyle',opt.LINESTYLE{i},...
                        'linewidth',opt.LINEWIDTH);
                end
            otherwise
                if(opt.ABSOLUTE)
                    rh=plot(opt.AXIS,...
                        z6(i)+b(i)/86400+...
                        (0:delta(i)/86400:delta(i)/86400*(npts(i)-1)).',...
                        data(i).dep,...
                        'color',opt.CMAP(i,:),...
                        'linestyle',opt.LINESTYLE{i},...
                        'linewidth',opt.LINEWIDTH);
                else
                    rh=plot(opt.AXIS,...
                        b(i)+(0:delta(i):delta(i)*(npts(i)-1)).',...
                        data(i).dep,...
                        'color',opt.CMAP(i,:),...
                        'linestyle',opt.LINESTYLE{i},...
                        'linewidth',opt.LINEWIDTH);
                end
        end
    end
    
    % add 1st handle to setup for markers
    userdata.markers.records(i)=rh(1);
    
    % set userdata to everything but the data (cleared)
    if(goodfiles(i))
        if(isempty(data(i).dep)); data(i).empty=true; end
        data(i).dep=[];
        data(i).ind=[];
    end
    
    % save info to handle
    nrh=numel(rh);
    for ridx=1:nrh
        if(nrh>1)
            ncmpstr=[' (CMP#' num2str(ridx) ')'];
        else
            ncmpstr=[];
        end
        data(i).index=[i ridx];
        data(i).yrng=yrng(i,:);
        set(rh(ridx),'userdata',data(i),'tag','record',...
            'displayname',[displayname{i} ncmpstr]);
    end
end
hold(opt.AXIS,'off');

% extras
box(opt.AXIS,'on');
grid(opt.AXIS,'on');
axis(opt.AXIS,'tight');

% add markers to axis userdata
userdata.markers.names=marknames;
userdata.markers.times=marktimes;
userdata.function='plot0';
set(opt.AXIS,'userdata',userdata);

% draw markers
if(opt.MARKERS)
    drawmarkers(opt.AXIS,varargin{:});
end

% special yaxis tick labels (names)
if(~isempty(opt.NAMESONYAXIS) && any(opt.NAMESONYAXIS))
    set(opt.AXIS,'ytick',1:nrecs);
    switch opt.NAMESONYAXIS
        case {'kstnm' 'st' 'sta' 'station'}
            set(opt.AXIS,'yticklabel',kname(:,2));
        case {'kcmpnm' 'cmp' 'component'}
            set(opt.AXIS,'yticklabel',kname(:,4));
        case {'shortcmp' 'short'}
            kname=char(kname(:,4));
            kname=cellstr(kname(:,3));
            set(opt.AXIS,'yticklabel',kname);
        case {'stshort' 'stashort' 'stationshort'}
            tmp=char(kname(:,4));
            tmp=cellstr(tmp(:,3));
            kname=strcat(kname(:,2),'.',tmp);
            set(opt.AXIS,'yticklabel',kname);
        case {'stcmp' 'stacmp'}
            kname=strcat(kname(:,2),'.',kname(:,4));
            set(opt.AXIS,'yticklabel',kname);
        case {'kname'}
            set(opt.AXIS,'yticklabel',displayname);
        otherwise
            set(opt.AXIS,'yticklabel',{data.name}.');
    end
    
    % make sure we can see all names
    ylimits=ylim(opt.AXIS);
    if(ylimits(1)>1); ylimits(1)=1; end
    if(ylimits(2)<nrecs); ylimits(2)=nrecs; end
    ylim(opt.AXIS,ylimits);
end

% axis zooming
if(~isempty(opt.XLIM))
    %if(isempty(opt.YLIM)); axis autoy; end
    xlim(opt.AXIS,opt.XLIM);
end
if(~isempty(opt.YLIM))
    %if(isempty(opt.XLIM)); axis autox; end
    ylim(opt.AXIS,opt.YLIM);
end

% datetick
if(opt.ABSOLUTE)
    if(isempty(opt.DATEFORMAT))
        if(isempty(opt.XLIM))
            datetick(opt.AXIS,'x');
        else
            datetick(opt.AXIS,'x','keeplimits');
        end
    else
        if(isempty(opt.XLIM))
            datetick(opt.AXIS,'x',opt.DATEFORMAT);
        else
            datetick(opt.AXIS,'x',opt.DATEFORMAT,'keeplimits');
        end
    end
else
    if(~isempty(opt.DATEFORMAT))
        if(isempty(opt.XLIM))
            datetick(opt.AXIS,'x',opt.DATEFORMAT);
        else
            datetick(opt.AXIS,'x',opt.DATEFORMAT,'keeplimits');
        end
    end
end

% label
if(~isempty(opt.TITLE) && isnumeric(opt.TITLE) && opt.TITLE==1)
    opt.TITLE=[num2str(numel(goodfiles)) ...
        '/' num2str(nrecs) ' Records'];
end
if(~isempty(opt.XLABEL) && isnumeric(opt.XLABEL) && opt.XLABEL==1)
    if(opt.ABSOLUTE)
        xlimits=get(opt.AXIS,'xlim');
        opt.XLABEL=joinwords(cellstr(datestr(unique(fix(xlimits)))),...
            '   to   ');
    else
        opt.XLABEL='Time (sec)';
    end
end
if(~isempty(opt.YLABEL) && isnumeric(opt.YLABEL) && opt.YLABEL==1)
    opt.YLABEL='Record';
end
title(opt.AXIS,opt.TITLE,'color',opt.FGCOLOR,'fontsize',opt.FONTSIZE);
xlabel(opt.AXIS,opt.XLABEL,'color',opt.FGCOLOR,'fontsize',opt.FONTSIZE);
ylabel(opt.AXIS,opt.YLABEL,'color',opt.FGCOLOR,'fontsize',opt.FONTSIZE);

% output axes if wanted
if(nargout); varargout{1}=opt.AXIS; end

end
