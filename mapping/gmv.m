function [varargout]=gmv(data,varargin)
%GMV    Ground Motion Visualization via an evolving station map
%
%    Usage:    gmv(data)
%              gmv(data,'prop1',val1,'prop2',val2,...)
%              mov=gmv(...)
%
%    Description:
%     GMV(DATA) creates a ground motion visualization (gmv) using a station
%     map where the station elements evolve with time based on the ground
%     motion they are experiencing.  For vertical motion a colored marker
%     changes from red (positive) to blue (negative).  For horizontal
%     motion a radial black line indicates the direction and magnitude at
%     the corresponding location.  DATA is a SEIZMO data struct where all
%     of the data shares the same timing for each sample.  You can use
%     INTERPOLATE+SYNCHRONIZE to assure this (see Examples section below).
%
%     GMV(DATA,'PROP1',VAL1,'PROP2',VAL2,...) alters certain aspects of the
%     gmv by adjusting property values.  Available properties are:
%      'mode'    - 'vert', 'horz', 'both' (default is 'both')
%      'vcmap'   - colormap for vertical ground motions (def is 'blue2red')
%      'vclip'   - limits in vert ground motion (2 std)
%      'hclip'   - limits in horz ground motion (2 std)
%      'hscale'  - horz ground motion that gives a 1deg bar (1 std)
%      'vstd'    - number of vert std dev for default vclip (2)
%      'hstd'    - number of horz std dev for default hclip (2)
%      'hstdsc'  - number of horz std dev for default hscale (1)
%      'hlinew'  - linewidth for horizontal motion line (def is 2)
%      'rsta'    - ref stn idx (in DATA) or 'name' (default is 1)
%      'xfunc'   - function to call to alter plot just before frame capture
%      'mmapopt' - options passed to mmap (put them in a cell array)
%      'p0opt'   - options passed to plot0 (put them in a cell array)
%      'delay'   - delay (in seconds) between frames (default is 0)
%
%     MOV=GMV(...) returns the movie.  This is useful for repeated
%     visualization (using MOVIE) or saving to an .avi file (using
%     MOVIE2AVI).
%
%    Notes:
%     - Based on the pretty movies by Chuck Ammon & Bob Woodward
%     - see the code for help with the 'xfunc' option
%
%    Examples:
%     % GMV requires the dataset to be synchronized.  Here we sample a
%     % dataset at 1Hz for the first hour after an event:
%     mov=gmv(interpolate(synchronize(data,'o',[],'io'),1,[],0,3600,0));
%
%     % Make a movie and save to an .avi file:
%     mov=gmv(data);
%     movie2avi(mov,'gmv.avi');
%     unixcompressavi('gmv.avi'); % for Unix/Linux w/ Mencoder installed
%     unix('mplayer gmv.avi');
%
%    See also: MAPSTATIONS, MMAP, MOVIE, MOVIE2AVI, UNIXCOMPRESSAVI

%     Version History:
%        Apr. 16, 2011 - initial version
%        Apr.  3, 2012 - use seizmocheck
%        Aug. 31, 2012 - handle rotate breaking if no horizontals, allow
%                        setting labels on record axis, minor doc update,
%                        handle absolute time plot
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 31, 2012 at 14:05 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));
if(nargin>1 && ~mod(nargin,2))
    error('seizmo:gmv:propsMustBePaired',...
        '1 or more properties were not paired with a value!');
elseif(nargin>1 && ~iscellstr(varargin(1:2:end)))
    error('seizmo:gmv:badProps',...
        '1 or more properties were not given as strings!');
end

% check data structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    data=checkheader(data,...
        'MULCMP_DEP','ERROR',...
        'NONTIME_IFTYPE','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'MULTIPLE_DELTA','ERROR',...
        'MULTIPLE_NPTS','ERROR',...
        'NONINTEGER_REFTIME','ERROR',...
        'UNSET_REFTIME','ERROR',...
        'OUTOFRANGE_REFTIME','ERROR',...
        'UNSET_ST_LATLON','ERROR');
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% make gmv
try
    % verbosity
    verbose=seizmoverbose;
    
    % get necessary header info
    [npts,delta,kstnm]=getheader(data,'npts','delta','kstnm');
    npts=npts(1);
    delta=delta(1);
    
    % parse properties
    prop=parse_gmv_props(varargin{:});
    if(isa(prop.vcmap,'function_handle'))
        prop.vcmap=prop.vcmap(64);
    end
    
    % get indices for reference station records
    if(ischar(prop.rsta))
        prop.rsta=find(strcmpi(prop.rsta,kstnm));
    end
    
    % vert & horz components
    switch prop.mode
        case 'vert'
            vdata=data(vertcmp(data));
            hdata=[];
            if(isempty(vdata))
                error('seizmo:gmv:noVert',...
                    'No vertical data found to visualize!');
            end
        case 'horz'
            vdata=[];
            try
                hdata=rotate(data,'to',0);
            catch
                error('seizmo:gmv:noHorz',...
                    'No horizontal pairs found to visualize!');
            end
        case 'both'
            vdata=data(vertcmp(data));
            try
                hdata=rotate(data,'to',0);
            catch
                hdata=[];
            end
            if(isempty(vdata) && isempty(hdata))
                error('seizmo:gmv:noHorz',...
                    'No vertical/horizontal data found to visualize!');
            end
            if(isempty(vdata)); prop.mode='horz'; end
            if(isempty(hdata)); prop.mode='vert'; end
    end
    
    % auto amplitude scaling
    seizmoverbose(false);
    switch prop.mode
        case {'both' 'vert'}
            if(isempty(prop.vclip))
                prop.vclip=mean(getvaluefun(vdata,@std));
                prop.vclip=prop.vstd*[-1 1]*prop.vclip;
            end
    end
    switch prop.mode
        case {'both' 'horz'}
            if(isempty(prop.hclip))
                prop.hclip=prop.hstd*mean(getvaluefun(hdata,@std));
            end
            if(isempty(prop.hscale))
                prop.hscale=prop.hstdsc*mean(getvaluefun(hdata,@std));
            end
    end
    seizmoverbose(verbose);
    
    % detail message
    if(verbose)
        disp('Making Ground Motion Visualization');
        print_time_left(0,npts);
    end
    
    % create gmv plot
    [fh,mapax,recax,xdata]=plot_gmv(data,vdata,hdata,prop);
    
    % save frame
    if(nargout); varargout{1}(1)=getframe(fh); end
    
    % variables for loop
    switch prop.mode
        case {'both' 'vert'}
            vh=findobj(mapax,'tag','stations');
    end
    switch prop.mode
        case {'both' 'horz'}
            [la1,lo1]=getheader(hdata(1:2:end),'stla','stlo');
    end
    
    % update loop
    for i=2:npts
        % forced delay
        drawnow;
        pause(prop.delay);
        
        % quit if closed
        if(~ishandle(fh)); return; end
        
        % update map
        switch prop.mode
            case {'vert' 'both'}
                z=singlepnt(vdata,i);
                c=z2c(z,prop.vcmap,prop.vclip);
                set(vh,'cdata',c);
                drawnow;
        end
        switch prop.mode
            case {'horz' 'both'}
                % get particle amplitudes & azimuths
                % plot as great circles (works in polar regions)
                delete(findobj(mapax,'tag','horzgm'));
                n=singlepnt(hdata(1:2:end),i);
                e=singlepnt(hdata(2:2:end),i);
                az=90-atan2(n,e)*180/pi;
                gc=min(prop.hclip,sqrt(n.^2+e.^2))/prop.hscale;
                [la2,lo2]=sphericalfwd(la1,lo1,gc,az);
                [lat,lon]=gcarc2latlon(la1,lo1,la2,lo2,5);
                lon=unwrap(lon*pi/180,[],2)*180/pi; % avoid wraparound
                m_line(lon',lat','linewi',prop.hlinew,...
                    'color','k','tag','horzgm');
        end
        switch prop.mode
            case 'both'
                movekids(vh,'front');
        end
        movekids(findobj(mapax,'tag','rstaloc'),'front');
        
        % update record plot
        set(findobj(recax,'tag','currenttime'),'xdata',xdata(i)*[1;1]);
        
        % save frame
        if(nargout); varargout{1}(i)=getframe(fh); end
        
        % detail message
        if(verbose); print_time_left(i,npts); end
    end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror);
end

end


function [prop]=parse_gmv_props(varargin)
% defaults
%      'mode'    - 'vert', 'horz', 'both' (default is 'both')
%      'vcmap'   - colormap for vertical ground motions (def is 'blue2red')
%      'vclip'   - limits in vert ground motion (def set to 2 std)
%      'hclip'   - limit in horz ground motion (def set to 2 std)
%      'hscale'  - horz ground motion that gives a 1deg bar (1 std)
%      'vstd'    - number of vert standard deviations for vclip (2)
%      'hstd'    - number of horz standard deviations for hclip (2)
%      'hstdsc'  - number of horz standard deviations for hscale (1)
%      'hlinew'  - linewidth for horizontal motion line (def is 2)
%      'rsta'    - ref stn idx (in DATA) or 'name' (default is 1)
%      'xfunc'   - function to call to alter plot just before frame capture
%      'mmapopt' - options passed to mmap (put them in a cell array)
%      'p0opt'   - options passed to plot0 (put them in a cell array)
%      'delay'   - delay (in seconds) between frames (default is 0)
prop.mode='both';
prop.vcmap=@blue2red;
prop.vclip=[];
prop.hclip=[];
prop.hscale=[];
prop.vstd=2;
prop.hstd=2;
prop.hstdsc=1;
prop.hlinew=2;
prop.rsta=1;
prop.xfunc=@xfunc;
prop.mmapopt={};
prop.p0opt={};
prop.delay=0;

for i=1:2:nargin
    switch lower(varargin{i})
        case 'mode'
            prop.mode=varargin{i+1};
        case {'vcmap' 'vertcmap' 'vcolormap' 'vertcolormap'}
            prop.vcmap=varargin{i+1};
        case {'vclip' 'vertclip' 'vlim' 'vertlim' 'vlimit' 'vlimits'}
            prop.vclip=varargin{i+1};
        case {'hclip' 'horzclip' 'hlim' 'horzlim' 'hlimit' 'hlimits'}
            prop.hclip=varargin{i+1};
        case {'hscale' 'horzscale'}
            prop.hscale=varargin{i+1};
        case {'vstd' 'vstddev'}
            prop.vstd=varargin{i+1};
        case {'hstd' 'hstddev'}
            prop.hstd=varargin{i+1};
        case {'hstdsc' 'hstdscale' 'hstddevsc' 'hstddevscale'}
            prop.hstdsc=varargin{i+1};
        case {'hlinew' 'hlinewidth'}
            prop.hlinew=varargin{i+1};
        case {'rsta' 'rst' 'refst' 'refsta' 'ref' 'r' 'st' 'sta'}
            prop.rsta=varargin{i+1};
        case {'xfunc' 'x' 'xfun' 'xf' 'f' 'fu' 'fun' 'func'}
            prop.xfunc=varargin{i+1};
        case {'mmapopt' 'mo' 'mopt' 'mmo' 'mmopt' 'mmap'}
            prop.mmapopt=varargin{i+1};
        case {'p0opt' 'po' 'popt' 'p0'}
            prop.p0opt=varargin{i+1};
        case {'delay' 'd' 'del'}
            prop.delay=varargin{i+1};
        otherwise
            error('seizmo:gmv:badInput',...
                'Unknown property: %s',varargin{i});
    end
end

end


function [amp]=singlepnt(data,idx)
nrecs=numel(data);
amp=nan(nrecs,1);
for i=1:nrecs; amp(i)=data(i).dep(idx); end
end


function [fh,mapax,recax,xdata]=plot_gmv(data,vdata,hdata,prop)
% initialize
fh=figure('color','k');
mapax=subplot(4,1,1:3,'parent',fh);
recax=subplot(4,1,4,'parent',fh);

% mmap
switch prop.mode
    case {'vert' 'both'}
        mapax=mapstations(vdata,'axis',mapax,prop.mmapopt{:});
        z=singlepnt(vdata,1);
        c=z2c(z,prop.vcmap,prop.vclip);
        set(findobj(mapax,'tag','stations'),'cdata',c);
        drawnow;
    case 'horz'
        mapax=mmap('axis',mapax,prop.mmapopt{:});
end
switch prop.mode
    case {'horz' 'both'}
        % get particle amplitudes & azimuths
        % plot as great circles (works in polar regions)
        [la1,lo1]=getheader(hdata(1:2:end),'stla','stlo');
        n=singlepnt(hdata(1:2:end),1);
        e=singlepnt(hdata(2:2:end),1);
        az=90-atan2(n,e)*180/pi;
        gc=min(prop.hclip,sqrt(n.^2+e.^2))/prop.hscale;
        [la2,lo2]=sphericalfwd(la1,lo1,gc,az);
        [lat,lon]=gcarc2latlon(la1,lo1,la2,lo2,5);
        lon=unwrap(lon*pi/180,[],2)*180/pi; % avoid wraparound streak
        m_line(lon',lat','linewi',prop.hlinew,'color','k','tag','horzgm');
end
switch prop.mode
    case 'both'
        movekids(findobj(mapax,'tag','stations'),'front');
end

% highlight reference station location
hold(mapax,'on');
[lat,lon]=getheader(data(prop.rsta),'stla','stlo');
m_scatter(mapax,lon,lat,[],'y','tag','rstaloc');
hold(mapax,'off');

% record plot
recax=plot0(data(prop.rsta),'axis',recax,'namesonyaxis','stashort',...
    'title',[],'ylabel',[],'xlabel','Time Since Earthquake (seconds)',...
    prop.p0opt{:});
hold(recax,'on');

% add current time bar
ylimits=ylim(recax);
kids=get(recax,'children');
xdata=get(kids(1),'xdata');
plot(recax,xdata(1)*[1;1],ylimits.','r','linewidth',2,'tag','currenttime');
hold(recax,'off');

% modify using external function
prop.xfunc(fh,mapax,recax,prop);
end


function []=xfunc(varargin)
% does nothing
end

