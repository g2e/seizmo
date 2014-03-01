function [varargout]=plot2(varargin)
%PLOT2    Overlay plot of SEIZMO data records
%
%    Usage:    plot2(data)
%              plot2(data1,data2) or plot2(data1,...,dataN)
%              plot2(...,'option',value,...)
%              ax=plot2(...)
%
%    Description:
%     PLOT2(DATA) draws all non-xyz records in SEIZMO struct DATA over one
%     another in the same plot in a new figure.  Each record is plotted as
%     a distinct color from the HSV colormap.  Spectral records are
%     converted to the time domain prior to plotting.
%
%     PLOT2(DATA1,DATA2) OR PLOT2(DATA1,...,DATAN) will draw records with
%     the same index in their dataset in the same subplot.  So DATA1(3) and
%     DATA2(3) will be plotted together in the 3rd subplot.  All datasets
%     must have the same number of records or be scalar.  This is extremely
%     useful for plotting real vs synthetic data.
%
%     PLOT2(...,'OPTION',VALUE,...) sets certain plotting options to do
%     simple manipulation of the plots.  Available options are:
%      FGCOLOR    -- foreground color (axes, text, labels)
%      BGCOLOR    -- background color (does not set figure color)
%      AXIS       -- axes to plot in (must be multiple if multidataset)
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
%
%     AX=PLOT2(...) returns the handles for all the axes drawn in.  This is
%     useful for more detailed plot manipulation.
%
%    Notes:
%
%    Examples:
%     % to overlay the first 4 records in a dataset:
%     plot2(data(1:4))
%
%     % same but in UTC time
%     plot2(data(1:4),'utc',true)
%
%     % plot first 4 records against the next 4
%     plot2(data(1:4),data(5:8))
%
%     %o overlay the 5th and 8th records from 0 to 300 seconds
%     plot2(data([5 8]),'xlim',[0 300]);
%
%    See also: PLOT0, PLOT1, RECORDSECTION

%     Version History:
%        Aug. 14, 2010 - rewrite
%        Sep. 16, 2010 - tagged record handles as 'record'
%        Feb. 10, 2011 - clear todo list (use datenum('datetime') to do
%                        absolute time limits with strings)
%        Apr. 16, 2011 - allow empty title/xlabel/ylabel, allow datetick
%                        for non-absolute
%        Apr. 19, 2011 - userdata for each record contains record metadata
%        Feb.  6, 2012 - better getheader usage, set displayname to kname
%                        string, fix legend coloring
%        Aug. 30, 2012 - allow [] input for labeling
%        Mar.  1, 2014 - maketime fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 23:00 GMT

% todo:
% - bug: kname output for title in multi-dataset case only for last dataset

% check nargin
error(nargchk(1,inf,nargin));

% find seizmo structures
issz=false(nargin,1);
for i=1:nargin
    issz(i)=isseizmo(varargin{i},'dep');
end
nd=sum(issz);
if(~nd)
    error('seizmo:plot2:badInput',...
        'No datasets provided!');
end

% separate datasets from option/value pairs
data=varargin(issz);
varargin(issz)=[];

% require all datasets to have equal number of records
if(~isequalnumelorscalar(data{:}))
    error('seizmo:plot2:badInput',...
        'Datasets must have the same number of records!');
end

% flatten datasets and expand single record datasets
[data{:}]=mat2vec(data{:});
[data{:}]=expandscalars(data{:});
nrecs=numel(data{1});

% default/parse options
opt=parse_seizmo_plot_options(varargin{:});

% line coloring
n=nrecs; if(nd>1); n=nd; end
if(ischar(opt.CMAP) || iscellstr(opt.CMAP))
    % list of color names or a colormap function
    try
        % attempt color name to rgb conversion first
        opt.CMAP=name2rgb(opt.CMAP);
        opt.CMAP=repmat(opt.CMAP,ceil(n/size(opt.CMAP,1)),1);
    catch
        % guess its a colormap function then
        opt.CMAP=str2func(opt.CMAP);
        opt.CMAP=opt.CMAP(n);
    end
else
    % numeric colormap array
    opt.CMAP=repmat(opt.CMAP,ceil(n/size(opt.CMAP,1)),1);
end

% line style
opt.LINESTYLE=cellstr(opt.LINESTYLE);
opt.LINESTYLE=opt.LINESTYLE(:);
opt.LINESTYLE=repmat(opt.LINESTYLE,ceil(n/size(opt.LINESTYLE,1)),1);

% check axes handle
if(nd>1)
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
    
    % necessary header info
    [leven,displayname]=deal(cell(nrecs,nd));
    [b,npts,delta,z]=deal(nan(nrecs,nd));
    goodfiles=false(nrecs,nd);
    for i=1:nd
        % check filetype
        iftype=getheader(data{i},'iftype id');
        time=strcmpi(iftype,'itime') | strcmpi(iftype,'ixy');
        spec=strcmpi(iftype,'irlim') | strcmpi(iftype,'iamph');
        goodfiles(:,i)=time | spec;
        
        % convert spectral to timeseries
        if(sum(spec)); data{i}(spec)=idft(data{i}(spec)); end
        
        % header info
        [b(:,i),npts(:,i),delta(:,i),kname,z6,leven(:,i)]=...
            getheader(data{i},'b','npts','delta','kname','z6','leven lgc');
        z(:,i)=datenum(cell2mat(z6));
        
        % names for legend
        displayname(:,i)=strcat(kname(:,1),'.',kname(:,2),...
            '.',kname(:,3),'.',kname(:,4));
    end
    
    % loop over each subplot
    for i=1:nrecs
        % adjust current axis
        cla(opt.AXIS(i),'reset');
        set(opt.AXIS(i),'ydir',opt.YDIR,'xdir',opt.XDIR,...
            'xscale',opt.XSCALE,'yscale',opt.YSCALE,...
            'fontsize',opt.FONTSIZE,'color',opt.BGCOLOR,...
            'xcolor',opt.FGCOLOR,'ycolor',opt.FGCOLOR);
        
        % loop over datasets
        hold(opt.AXIS(i),'on');
        for j=1:nd
            if(~goodfiles(i,j)); continue; end
            switch leven{i}
                case 'false'
                    if(opt.ABSOLUTE)
                        rh=plot(opt.AXIS(i),...
                            z(i,j)+data{j}(i).ind/86400,data{j}(i).dep,...
                            'color',opt.CMAP(j,:),...
                            'linestyle',opt.LINESTYLE{j},...
                            'linewidth',opt.LINEWIDTH);
                    else
                        rh=plot(opt.AXIS(i),data{j}(i).ind,data{j}(i).dep,...
                            'color',opt.CMAP(j,:),...
                            'linestyle',opt.LINESTYLE{j},...
                            'linewidth',opt.LINEWIDTH);
                    end
                otherwise
                    if(opt.ABSOLUTE)
                        rh=plot(opt.AXIS(i),z(i,j)+b(i,j)/86400+...
                            (0:delta(i,j)/86400:delta(i,j)/86400*...
                            (npts(i,j)-1)).',data{j}(i).dep,...
                            'color',opt.CMAP(j,:),...
                            'linestyle',opt.LINESTYLE{j},...
                            'linewidth',opt.LINEWIDTH);
                    else
                        rh=plot(opt.AXIS(i),b(i,j)+...
                            (0:delta(i,j):delta(i,j)*(npts(i,j)-1)).',...
                            data{j}(i).dep,...
                            'color',opt.CMAP(j,:),...
                            'linestyle',opt.LINESTYLE{j},...
                            'linewidth',opt.LINEWIDTH);
                    end
            end
            
            % set userdata to everything but data (cleared)
            data{j}(i).dep=[];
            data{j}(i).ind=[];
            nrh=numel(rh);
            for ridx=1:nrh
                if(nrh>1)
                    ncmpstr=[' (' num2str(ridx) ')'];
                else
                    ncmpstr=[];
                end
                data{j}(i).index=[j i ridx];
                set(rh(ridx),'userdata',data{j}(i),...
                    'displayname',[displayname{i,j} ncmpstr]);
            end
        end
        hold(opt.AXIS(i),'off');
        
        % tag records
        rh=get(opt.AXIS(i),'children');
        set(rh,'tag','record');
        
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
                    datetick(opt.AXIS(i),'x','keeplimits','keepticks');
                else
                    datetick(opt.AXIS(i),'x','keeplimits','keepticks');
                end
            else
                if(isempty(opt.XLIM))
                    datetick(opt.AXIS(i),'x',opt.DATEFORMAT,...
                        'keeplimits','keepticks');
                else
                    datetick(opt.AXIS(i),'x',opt.DATEFORMAT,...
                        'keeplimits','keepticks');
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
        % - FIXME: KNAME belongs to last dataset!
        if(~isempty(opt.TITLE) && isnumeric(opt.TITLE))
            switch opt.TITLE
                case 1 % x/y records
                    p2title=[num2str(sum(goodfiles(i,:))) ...
                        '/' num2str(nd) ' Records'];
                case 2 % kstnm
                    p2title=kname(i,2);
                case 3 % kcmpnm
                    p2title=kname(i,4);
                case 4 % shortcmp
                    p2title=kname{i,4}(3);
                case 5 % stashort
                    p2title=strcat(kname(i,2),'.',kname{i,4}(3));
                case 6 % stcmp
                    p2title=strcat(kname(i,2),'.',kname(i,4));
                case 7 % kname
                    p2title=texlabel(strcat(kname(i,1),'.',kname(i,2),...
                        '.',kname(i,3),'.',kname(i,4)),'literal');
                otherwise
                    p2title=['RECORD ' num2str(i)];
            end
        else
            p2title=opt.TITLE;
        end
        if(~isempty(opt.XLABEL) && isnumeric(opt.XLABEL) && opt.XLABEL==1)
            if(opt.ABSOLUTE)
                xlimits=get(opt.AXIS(i),'xlim');
                p2xlabel=...
                    joinwords(cellstr(datestr(unique(fix(xlimits)))),...
                    '  to  ');
            else
                p2xlabel='Time (sec)';
            end
        else
            p2xlabel=opt.XLABEL;
        end
        if(~isempty(opt.YLABEL) && isnumeric(opt.YLABEL) && opt.YLABEL==1)
            p2ylabel='Amplitude';
        else
            p2ylabel=opt.YLABEL;
        end
        title(opt.AXIS(i),p2title,...
            'color',opt.FGCOLOR,'fontsize',opt.FONTSIZE);
        xlabel(opt.AXIS(i),p2xlabel,...
            'color',opt.FGCOLOR,'fontsize',opt.FONTSIZE);
        ylabel(opt.AXIS(i),p2ylabel,...
            'color',opt.FGCOLOR,'fontsize',opt.FONTSIZE);
    end
else
    % all in one plot
    if(isempty(opt.AXIS) || ~isscalar(opt.AXIS) || ~isreal(opt.AXIS) ...
            || ~ishandle(opt.AXIS) || ~strcmp('axes',get(opt.AXIS,'type')))
        % new figure
        figure('color',opt.BGCOLOR,'defaulttextcolor',opt.FGCOLOR,...
            'defaultaxesxcolor',opt.FGCOLOR); % defaults for legend
        opt.AXIS=gca;
    else
        cla(opt.AXIS,'reset');
    end
    
    % adjust current axis
    set(opt.AXIS,'ydir',opt.YDIR,'xdir',opt.XDIR,...
        'xscale',opt.XSCALE,'yscale',opt.YSCALE,...
        'fontsize',opt.FONTSIZE,'color',opt.BGCOLOR,...
        'xcolor',opt.FGCOLOR,'ycolor',opt.FGCOLOR);
    
    % decell dataset
    data=data{1};
    
    % check filetype
    iftype=getheader(data,'iftype id');
    time=strcmpi(iftype,'itime') | strcmpi(iftype,'ixy');
    spec=strcmpi(iftype,'irlim') | strcmpi(iftype,'iamph');
    goodfiles=find(time | spec)';
    
    % convert spectral to timeseries
    if(sum(spec)); data(spec)=idft(data(spec)); end
    
    % necessary header info
    [b,npts,delta,z6,kname,leven]=getheader(data,...
        'b','npts','delta','z6','kname','leven lgc');
    z6=datenum(cell2mat(z6));
    
    % names for legend
    displayname=strcat(kname(:,1),'.',kname(:,2),...
        '.',kname(:,3),'.',kname(:,4));
    
    % loop through every record
    hold(opt.AXIS,'on');
    for i=goodfiles
        switch leven{i}
            case 'false'
                if(opt.ABSOLUTE)
                    rh=plot(opt.AXIS,z6(i)+data(i).ind/86400,...
                        data(i).dep,...
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
                    rh=plot(opt.AXIS,z6(i)+b(i)/86400+...
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
        opt.YLABEL='Amplitude';
    end
    title(opt.AXIS,opt.TITLE,...
        'color',opt.FGCOLOR,'fontsize',opt.FONTSIZE);
    xlabel(opt.AXIS,opt.XLABEL,...
        'color',opt.FGCOLOR,'fontsize',opt.FONTSIZE);
    ylabel(opt.AXIS,opt.YLABEL,...
        'color',opt.FGCOLOR,'fontsize',opt.FONTSIZE);
end

% output axes if wanted
if(nargout); varargout{1}=opt.AXIS; end

end
