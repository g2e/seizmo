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
%      FONTWEIGHT   -- 'light', 'normal', 'demi' or 'bold'
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
%        Apr. 19, 2011 - userdata for each record contains record metadata
%        Jan. 12, 2012 - minor improvement to normstyle handling
%        Jan. 25, 2012 - norm2yaxis bugfix
%        Feb.  6, 2012 - better getheader usage, set displayname to kname
%                        string, fix legend coloring, handle undefined
%                        yfield (errors/warns), fix ncmp>1 brokeness
%        Mar.  6, 2012 - don't use depmin/depmax from header
%        Mar. 13, 2012 - make sure dep* is updated for extract_plot_data
%        Aug. 30, 2012 - allow [] input for labeling
%        Mar.  1, 2014 - maketime fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 23:00 GMT

% todo:
% - bug: markers only based on 1st cmp data values

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
iftype=getheader(data,'iftype id');
time=strcmpi(iftype,'itime') | strcmpi(iftype,'ixy');
spec=strcmpi(iftype,'irlim') | strcmpi(iftype,'iamph');
goodfiles=(time | spec)';

% convert spectral to timeseries
if(sum(spec)); data(spec)=idft(data(spec)); end

% fix dep*
oldcheckheader_state=checkheader_state(true);
data=checkheader(data,'all','ignore','old_dep_stats','fix');
checkheader_state(oldcheckheader_state);

% header info
[b,npts,delta,z6,kname,leven,yfield]=getheader(data,'b',...
    'npts','delta','z6','kname','leven lgc',opt.YFIELD);
z6=datenum(cell2mat(z6));

% check for yfield=nan
bad=isnan(yfield);
if(all(bad))
    error('seizmo:recordsection:badInput',...
        'YFIELD (%s) for all records is undefined!',opt.YFIELD);
elseif(any(bad))
    warning('seizmo:recordsection:badInput',...
        ['Skipped plotting of record(s) with YFIELD (%s) undefined:\n' ...
        sprintf('%d ',find(bad))],opt.YFIELD);
    goodfiles=goodfiles & ~bad;
end
goodfiles=find(goodfiles);

% names for legend
displayname=strcat(kname(:,1),'.',kname(:,2),...
    '.',kname(:,3),'.',kname(:,4));

% get markers info
[marknames,marktimes]=getmarkers(data);

% convert markers to absolute time if used
if(opt.ABSOLUTE)
    marktimes=marktimes/86400+z6(:,ones(1,size(marktimes,2)));
end

% get data limits
depmin=nan(nrecs,1); depmax=nan(nrecs,1);
switch lower(opt.AMPSCALE)
    case 'linear'
        for i=goodfiles
            depmin(i)=abs(min(data(i).dep(:)));
            depmax(i)=abs(max(data(i).dep(:)));
        end
    case 'log'
        for i=goodfiles
            if(any(data(i).dep>0))
                depmin(i)=min(log10(data(i).dep(data(i).dep(:)>0)));
                depmax(i)=max(log10(data(i).dep(data(i).dep(:)>0)));
            end
        end
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
    scale=opt.NORMMAX;
end
switch opt.NORMSTYLE
    case {1 'i' 'single' 'individually' 'individual' 'one' 'separately'}
        switch lower(opt.AMPSCALE)
            case 'linear'
                ampmax=max(depmin,depmax);
                ampmax(ampmax==0)=1; % avoid divide-by-zero below
                for i=goodfiles
                    data(i).dep=yfield(i)+data(i).dep/ampmax(i)*scale;
                end
            case 'log'
                logrng=depmax-depmin;
                logrng(logrng==0)=1; % avoid divide-by-zero below
                for i=goodfiles
                    data(i).dep=yfield(i)+(2*((log10(data(i).dep)...
                        -depmin(i))/logrng(i))-1)*scale;
                end
        end
    case {0 'g' 'a' 'group' 'together' 'all'}
        switch lower(opt.AMPSCALE)
            case 'linear'
                ampmax=max([depmin; depmax]);
                ampmax(ampmax==0)=1; % avoid divide-by-zero below
                for i=goodfiles
                    data(i).dep=yfield(i)+data(i).dep/ampmax*scale;
                end
            case 'log'
                logmin=min(depmin);
                logmax=max(depmax);
                logrng=logmax-logmin;
                logrng(logrng==0)=1; % avoid divide-by-zero below
                for i=goodfiles
                    data(i).dep=yfield(i)+...
                        (2*((log10(data(i).dep)-logmin)/logrng)-1)*scale;
                end
        end
    otherwise
        error('seizmo:recordsection:badInput',...
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

% loop through every record
hold(opt.AXIS,'on');
for i=goodfiles
    switch leven{i}
        case 'false'
            if(opt.ABSOLUTE)
                rh=plot(opt.AXIS,z6(i)+data(i).ind/86400,data(i).dep,...
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
    
    % add 1st handle to setup for markers
    userdata.markers.records(i)=rh(1);
    
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

