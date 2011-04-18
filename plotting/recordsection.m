function [varargout]=recordsection(data,varargin)
%RECORDSECTION    Plots SEIZMO data records in a record section
%
%    Usage:    recordsection(data)
%              recordsection(...,'option',value,...)
%              ax=recordsection(...)
%
%    Description:
%     RECORDSECTION(DATA) draws all non-xyz records in SEIZMO struct DATA
%     spaced out by their 'gcarc' header field values (degree distance
%     from event location).  The records are normalized as a group with a
%     maximum amplitude range corresponding to a third of the y axis range.
%     Each record is drawn as a distinct color from the HSV colormap.
%     Spectral records are converted to the time domain prior to plotting.
%
%     RECORDSECTION(...,'OPTION',VALUE,...) sets certain plotting options
%     to do simple manipulation of the plots.  Available options are:
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
%      XDIR         -- 'normal' or 'reverse'
%      YDIR         -- 'normal' or 'reverse'
%      FONTSIZE     -- size of fonts in the axes
%      YFIELD       -- header field for y-axis positioning of records
%      XSCALE       -- 'linear' or 'log'
%      YSCALE       -- 'linear' or 'log'
%      AMPSCALE     -- 'linear' or 'log'
%      MARKERS      -- true/false where true draws the markers
%
%     AX=RECORDSECTION(...) returns the handles for all the axes drawn in.
%     This is useful for more detailed plot manipulation.
%
%    Notes:
%
%    Examples:
%     % make a azimuthal record section
%     recordsection(data,'yfield',az);
%
%     % put the record section in absolute time
%     recordsection(data,'utc',true);
%
%    See also: PLOT0, PLOT1, PLOT2

%     Version History:
%        Aug. 14, 2010 - rewrite
%        Sep. 14, 2010 - added marker support
%        Feb.  6, 2011 - fixed new figure/axes call
%        Feb. 10, 2011 - cleared todo list (2 axes is too problematic)
%        Apr. 16, 2011 - allow empty title/xlabel/ylabel, allow datetick
%                        for non-absolute, single record bugfix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 16, 2011 at 23:00 GMT

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

% check filetype
iftype=getenumid(data,'iftype');
time=strcmpi(iftype,'itime') | strcmpi(iftype,'ixy');
spec=strcmpi(iftype,'irlim') | strcmpi(iftype,'iamph');
goodfiles=find(time | spec)';

% convert spectral to timeseries
if(sum(spec)); data(spec)=idft(data(spec)); end

% header info
leven=getlgc(data,'leven');
[b,npts,delta,depmin,depmax,z6,kname,yfield]=getheader(data,...
    'b','npts','delta','depmin','depmax','z6','kname',opt.YFIELD);
depmin=abs(depmin);
depmax=abs(depmax);
z6=datenum(cell2mat(z6));

% get markers info
[marknames,marktimes]=getmarkers(data);

% convert markers to absolute time if used
if(opt.ABSOLUTE)
    marktimes=marktimes/86400+z6(:,ones(1,size(marktimes,2)));
end

% normalize
if(opt.NORM2YAXIS)
    scale=(max(yfield)-min(yfield))*opt.NORMMAX/2;
    if(scale==0) % single record fix
        if(~isempty(opt.YLIM))
            scale=(max(opt.YLIM)-min(opt.YLIM))*opt.NORMMAX/2;
        else
            scale=opt.NORMMAX/2;
        end
    end
else
    scale=P.NORMMAX;
end
switch opt.NORMSTYLE
    case {'single' 'individually' 'individual' 'one' 'separately'}
        switch lower(opt.AMPSCALE)
            case 'linear'
                ampmax=max(depmin,depmax);
                ampmax(ampmax==0)=1; % avoid NaNs
                for i=goodfiles
                    data(i).dep=yfield(i)+data(i).dep/ampmax(i)*scale;
                end
            case 'log'
                for i=goodfiles
                    logmin=min(log10(data(i).dep(data(i).dep>0)));
                    logmax=max(log10(data(i).dep(data(i).dep>0)));
                    logrng=logmax-logmin;
                    if(~logrng); logrng=1; end % avoid NaNs
                    data(i).dep=yfield(i)+...
                        (2*((log10(data(i).dep)-logmin)/logrng)-1)*scale;
                end
        end
    case {'group' 'together' 'all'}
        switch lower(opt.AMPSCALE)
            case 'linear'
                ampmax=max([depmin; depmax]);
                ampmax(ampmax==0)=1; % avoid NaNs
                for i=goodfiles
                    data(i).dep=yfield(i)+data(i).dep/ampmax*scale;
                end
            case 'log'
                logmin=nan(nrecs,1);
                logmax=nan(nrecs,1);
                for i=goodfiles
                    logmin(i)=min(log10(data(i).dep(data(i).dep>0)));
                    logmax(i)=max(log10(data(i).dep(data(i).dep>0)));
                end
                logmin=min(logmin);
                logmax=max(logmax);
                logrng=logmax-logmin;
                if(~logrng); logrng=1; end % avoid NaNs
                for i=goodfiles
                    data(i).dep=yfield(i)+...
                        (2*((log10(data(i).dep)-logmin)/logrng)-1)*scale;
                end
        end
end

% all in one plot
if(isempty(opt.AXIS) || ~isscalar(opt.AXIS) || ~isreal(opt.AXIS) ...
        || ~ishandle(opt.AXIS) || ~strcmp('axes',get(opt.AXIS,'type')))
    % new figure
    fh=figure('color',opt.BGCOLOR);
    opt.AXIS=axes('parent',fh);
else
    cla(opt.AXIS,'reset');
end

% adjust current axis
set(opt.AXIS,'ydir',opt.YDIR,'xdir',opt.XDIR,...
    'xscale',opt.XSCALE,'yscale',opt.YSCALE,...
    'fontsize',opt.FONTSIZE,'color',opt.BGCOLOR,...
    'xcolor',opt.FGCOLOR,'ycolor',opt.FGCOLOR);

% loop through every record
hold(opt.AXIS,'on');
for i=goodfiles
    switch leven{i}
        case 'false'
            if(opt.ABSOLUTE)
                plot(opt.AXIS,z6(i)+data(i).ind/86400,data(i).dep,...
                    'color',opt.CMAP(i,:),...
                    'linestyle',opt.LINESTYLE{i},...
                    'linewidth',opt.LINEWIDTH);
            else
                plot(opt.AXIS,data(i).ind,data(i).dep,...
                    'color',opt.CMAP(i,:),...
                    'linestyle',opt.LINESTYLE{i},...
                    'linewidth',opt.LINEWIDTH);
            end
        otherwise
            if(opt.ABSOLUTE)
                plot(opt.AXIS,...
                    z6(i)+(b(i)+(0:npts(i)-1)*delta(i)).'/86400,...
                    data(i).dep,...
                    'color',opt.CMAP(i,:),...
                    'linestyle',opt.LINESTYLE{i},...
                    'linewidth',opt.LINEWIDTH);
            else
                plot(opt.AXIS,(b(i)+(0:npts(i)-1)*delta(i)).',...
                    data(i).dep,...
                    'color',opt.CMAP(i,:),...
                    'linestyle',opt.LINESTYLE{i},...
                    'linewidth',opt.LINEWIDTH);
            end
    end
end
hold(opt.AXIS,'off');

% tag records
rh=get(opt.AXIS,'children');
set(rh,'tag','record');

% extras
box(opt.AXIS,'on');
grid(opt.AXIS,'on');
axis(opt.AXIS,'tight');

% add markers to axis userdata
userdata.markers.records=nan(nrecs,1);
userdata.markers.records(goodfiles)=flipud(rh);
userdata.markers.names=marknames;
userdata.markers.times=marktimes;
userdata.function='recordsection';
set(opt.AXIS,'userdata',userdata);

% draw markers
if(opt.MARKERS)
    drawmarkers(opt.AXIS,varargin{:});
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
if(isnumeric(opt.TITLE) && opt.TITLE==1)
    opt.TITLE=[num2str(numel(goodfiles)) ...
        '/' num2str(nrecs) ' Records'];
end
if(isnumeric(opt.XLABEL) && opt.XLABEL==1)
    if(opt.ABSOLUTE)
        xlimits=get(opt.AXIS,'xlim');
        opt.XLABEL=joinwords(cellstr(datestr(unique(fix(xlimits)))),...
            '   to   ');
    else
        opt.XLABEL='Time (sec)';
    end
end
if(isnumeric(opt.YLABEL) && opt.YLABEL==1)
    switch lower(opt.YFIELD)
        case 'gcarc'
            opt.YLABEL='Distance (degrees)';
        case 'dist'
            opt.YLABEL='Distance (km)';
        case 'az'
            opt.YLABEL='Azimuth (degrees)';
        case 'baz'
            opt.YLABEL='Backazimuth (degrees)';
        otherwise
            opt.YLABEL=opt.YFIELD;
    end
end
title(opt.AXIS,opt.TITLE,'color',opt.FGCOLOR,'fontsize',opt.FONTSIZE);
xlabel(opt.AXIS,opt.XLABEL,'color',opt.FGCOLOR,'fontsize',opt.FONTSIZE);
ylabel(opt.AXIS,opt.YLABEL,'color',opt.FGCOLOR,'fontsize',opt.FONTSIZE);

% output axes if wanted
if(nargout); varargout{1}=opt.AXIS; end

end

