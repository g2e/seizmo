function [info,xc,data0]=useralign(data,varargin)
%USERALIGN    Interactive alignment of a signal for SEIZMO records
%
%    Usage:    [info,xc,data0]=useralign(data)
%              [...]=useralign(data,'option1',value1,...,'optionN',valueN)
%
%    Description:
%     [INFO,XC,DATA0]=USERALIGN(DATA) presents an interactive set of menus
%     & plots to guide the aligning of records in SEIZMO struct DATA on a
%     particular signal.  The workflow includes interactively applying a
%     moveout, windowing, tapering, adjusting correlation options.  Once
%     the solving is complete the user is allowed to pick the onset of the
%     stack to get the arrival times as close as possible to "picked"
%     times.  Afterwards, a plot presents the results and a menu asks the
%     user to decide to accept those results or redo the processing.  INFO
%     is a struct containing several substructs providing details for the
%     options used in each subfunction (USERMOVEOUT, USERWINDOW, USERTAPER,
%     CORRELATE, TTSOLVE) as well as the final solution within
%     INFO.SOLUTION.  XC is the struct from CORRELATE reordered by TTSOLVE
%     (contains the correlogram peaks used to find the solution).  DATA0 is
%     the processed dataset (windowed, tapered & aligned).
%
%     [...]=USERALIGN(DATA,'OPTION1',VALUE1,...,'OPTIONN',VALUEN) passes
%     options to CUT, TAPER, CORRELATE & TTSOLVE for enhanced control of
%     the workflow.  Currently available are all options in TTSOLVE, 3
%     options in CORRELATE ('LAGS', 'NPEAKS', & 'ABSXC'), 'WINDOW' &
%     'TAPER.  'WINDOW' specifies the window passed to CUT and should be a
%     1x2 vector of [START END].  'TAPER' specifies the taper width passed
%     to TAPER.  See those functions for details about the options.  Note
%     that NPEAKS must be >=1!
%
%    Notes:
%
%    Examples:
%     % Say we have some recordings of a large teleseismic earthquake and
%     % that we have inserted the expected arrival times of Pdiff into the
%     % headers of these records.  Were did we get those expected arrival
%     % times?  The TauP routine TAUPTIME is used to raytrace through a 1D
%     % model such as PREM to get them (ARRIVALS2PICKS is probably the
%     % easiest method to get this done as it will call TAUPTIME and
%     % automatically add the arrival times to the header).  The first step
%     % is to "pre-align" on those times so that we can attempt to isolate
%     % the Pdiff waveforms from other seismic phases using windowing.  The
%     % Pdiff waveforms will have noticeable perturbations from this
%     % alignment because the 1D model does not account for the Earth's
%     % lateral heterogeneity in seismic properties.  To improve on this
%     % initial alignment, we can either go through and pick each signal
%     % onset manually (yuck!) or we can use waveform cross-correlation to
%     % solve for the best relative arrival times between all the recorded
%     % signals (there are other ways to align seismic phases but that is
%     % beyond the scope of this example).  USERALIGN provides us with an
%     % interface to do the later method.  It includes a set of menus and
%     % plots to aid us in both the signal isolation and alignment
%     % enhancing by solving for time perturbations utilizing multi-channel
%     % cross-correlation.  Great!  So this is pretty easy: pre-align using
%     % an initial guess such as expected arrivals from a 1D Earth model,
%     % isolate the signal using windowing and tapering, and finally
%     % realigning using cross-correlation:
%     data=timeshift(data,-findpicks(data,'Pdiff',1));
%     [info,xc,data0]=useralign(data);
%
%     % I usually like to give a SNR estimate of the waveforms in addition
%     % to the cross-correlation result to decide which waveforms are
%     % reliable:
%     data=timeshift(data,-findpicks(data,'Pdiff',1));
%     snr=quicksnr(data,[-100 -10],[-10 60]);
%     [info,xc,data0]=useralign(data,'snr',snr);
%
%    See also: TTSOLVE, CORRELATE, USERWINDOW, USERTAPER, USERRAISE,
%              USERMOVEOUT, USERALIGN_QUIET, MULTIBANDALIGN

%     Version History:
%        Mar. 16, 2010 - initial version (not really but whatever)
%        Mar. 17, 2010 - add xc, data0 to output, more crash options
%        Mar. 18, 2010 - output reordered xc (new TTSOLVE output), robust
%                        to menu closing
%        Mar. 22, 2010 - account for TTALIGN change, increase NPEAKS to 5
%        Mar. 24, 2010 - include stack picking to deal with dc-offset,
%                        user-driven cluster analysis for more info
%        Mar. 26, 2010 - drop cluster analysis (can be done separately),
%                        changed output (moved into info.solution)
%        Apr. 22, 2010 - replace crash with exit (but still crash)
%        Aug.  8, 2010 - doc update and slight code adjustments
%        Aug. 27, 2010 - update for new plotting functions, all in one
%                        figure
%        Sep. 15, 2010 - added code to better handle matched polarities,
%                        fixed new iteration bug (clear axes, not delete)
%        Oct.  2, 2010 - normalize amplitudes in plots
%        Nov. 13, 2010 - reverse order of aligned plot and stacking, update
%                        for userwindow arg change
%        Jan. 14, 2011 - add useralign_quiet to See also section
%        Jan. 17, 2011 - allow specifying window & taper, altered menus
%        Jan. 29, 2011 - fix window input bug
%        Mar. 16, 2012 - doc update
%        Oct. 17, 2012 - minor doc update (heterogeneity misspelling)
%        Jan. 11, 2013 - minor doc update (history update)
%        Jan. 30, 2013 - update for new correlate (no spacing option now
%                        but the lag option is now available)
%        Feb. 26, 2013 - bugfix: needs 'reltime' correlate option
%        Aug. 29, 2013 - bugfix: correlate settings switch cases were
%                        desynced from the menu options (thanks Cheng!)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 29, 2013 at 11:00 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% check data (dep)
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % check headers
    data=checkheader(data);
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% attempt align
try
    % basic checks on optional inputs
    if(mod(nargin-1,2))
        error('seizmo:useralign:badInput',...
            'Unpaired OPTION/VALUE!');
    elseif(~iscellstr(varargin(1:2:end)))
        error('seizmo:useralign:badInput',...
            'All OPTIONs must be specified with a string!');
    end
    
    % default correlate options
    info.correlate.lags=[];
    info.correlate.npeaks=5;
    info.correlate.absxc=true;
    info.userwindow.limits=[-inf inf];
    info.usertaper.width=[0 0];
    
    % parse out the correlate options
    keep=true(nargin-1,1);
    for i=1:2:nargin-1
        switch lower(varargin{i})
            case 'lags'
                if(numel(varargin{i+1})>2 || ~isnumeric(varargin{i+1}))
                    error('seizmo:useralign:badInput',...
                        'LAGS must be given as [MINLAG MAXLAG]!');
                end
                info.correlate.lags=varargin{i+1};
                keep(i:i+1)=false;
            case 'npeaks'
                if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}) ...
                        || varargin{i+1}~=fix(varargin{i+1}) ...
                        || varargin{i+1}<1)
                    error('seizmo:useralign:badInput',...
                        'NPEAKS must be an integer >=1 !');
                end
                info.correlate.npeaks=varargin{i+1};
                keep(i:i+1)=false;
            case 'absxc'
                if(~isscalar(varargin{i+1}) || ~islogical(varargin{i+1}))
                    error('seizmo:useralign:badInput',...
                        'ABSXC must be a logical value!');
                end
                info.correlate.absxc=varargin{i+1};
                keep(i:i+1)=false;
            case 'window'
                if(~isempty(varargin{i+1}) && (numel(varargin{i+1})~=2 ...
                        || ~isreal(varargin{i+1})))
                    error('seizmo:useralign_quiet:badInput',...
                        'WINDOW must be 1x2 vector as [START END]!');
                end
                info.userwindow.limits=varargin{i+1};
                keep(i:i+1)=false;
            case 'taper'
                if(~isempty(varargin{i+1}) && (~isreal(varargin{i+1}) ...
                        || numel(varargin{i+1})>2))
                    error('seizmo:useralign_quiet:badInput',...
                        'TAPER must be a scalar width like 0.05!');
                end
                info.usertaper.width=varargin{i+1};
                keep(i:i+1)=false;
        end
    end
    varargin=varargin(keep);
    
    % require the o field is set to the same value (within +/- 2 millisec)
    % - this essentially requires the recordings to be of a single quake
    o=getheader(data,'o');
    outc=cell2mat(getheader(data,'o utc'));
    if(any(abs(timediff(outc(1,:),outc))>0.002))
        error('seizmo:useralign:oFieldVaries',...
            'O field must correspond to one UTC time for all records!');
    end
    
    % copy data to working dataset (so user can go back)
    data0=data;
    
    % outer loop - only breaks free on user command
    happy_user=false;
    info.iteration=1;
    fh=figure('color','k');
    ax=makesubplots(2,3,1:5,...
        'parent',fh,'visible','off','color','k');
    while(~happy_user)
        % usermoveout
        set(ax(1),'visible','on');
        [data0,info.usermoveout,info.handles(1)]=usermoveout(data0,...
            'ax',ax(1),'normstyle','single');

        % userwindow
        set(ax(2),'visible','on');
        [data0,info.userwindow,info.handles(2)]=userwindow(data0,...
            info.userwindow.limits,[],[],'ax',ax(2),'normstyle','single');

        % usertaper
        set(ax(3),'visible','on');
        [data0,info.usertaper,info.handles(3)]=usertaper(data0,...
            info.usertaper.width,[],'ax',ax(3),'normstyle','single');

        % menu for correlate options
        while(1)
            % present current settings
            tmp1=num2str(info.correlate.npeaks);
            tmp2='POLARITIES ARE MATCHED';
            if(info.correlate.absxc); tmp2='POLARITIES ARE UNKNOWN'; end
            choice=menu('ALTER CORRELATE SETTINGS?',...
                ['CHANGE NUMBER OF PEAKS (' tmp1 ')'],...
                ['CHANGE POLARITY ASSUMPTION (' tmp2 ')'],...
                'NO - CORRELATE DATA');

            % proceed by user choice
            switch choice
                case 1 % npeaks
                    choice=menu('NUMBER OF PEAKS TO PICK',...
                        ['CURRENT (' tmp1 ')'],...
                        '1','3','5','7','9','CUSTOM');
                    switch choice
                        case 1 % CURRENT
                            % leave alone
                        case 2 % 1
                            info.correlate.npeaks=1;
                        case 3 % 3
                            info.correlate.npeaks=3;
                        case 4 % 5
                            info.correlate.npeaks=5;
                        case 5 % 7
                            info.correlate.npeaks=7;
                        case 6 % 9
                            info.correlate.npeaks=9;
                        case 7 % CUSTOM
                            tmp=inputdlg(...
                                ['Number of Peaks to Pick? [' tmp1 ']:'],...
                                'Custom Number of Peaks',1,{tmp1});
                            if(~isempty(tmp))
                                try
                                    tmp=str2double(tmp{:});
                                    if(isscalar(tmp) && isreal(tmp) ...
                                            && tmp==fix(tmp) && tmp>=1)
                                        info.correlate.npeaks=tmp;
                                    end
                                catch
                                    % do not change info.correlate.npeaks
                                end
                            end
                    end
                case 2 % absxc
                    choice=menu('DO THE POLARITIES ALL MATCH?',...
                        ['CURRENT (' tmp2 ')'],'YES','NO');
                    switch choice
                        case 1 % CURRENT
                            % leave alone
                        case 2 % looking at just positive peaks
                            info.correlate.absxc=false;
                        case 3 % looking at both pos/neg peaks
                            info.correlate.absxc=true;
                    end
                case 3 % correlate
                    break;
            end
        end

        % correlate
        if(info.correlate.absxc); absxc={'absxc'}; else absxc={}; end
        xc=correlate(data0,'mcxc','normxc','noauto','reltime',absxc{:},...
            info.correlate.lags,'peaks',{'npeaks',info.correlate.npeaks});

        % solve alignment
        if(info.correlate.absxc)
            [arr,err,pol,zmean,zstd,nc,info.ttsolve,xc]=ttsolve(xc,...
                varargin{:});
        else
            % force polarity solution
            [arr,err,pol,zmean,zstd,nc,info.ttsolve,xc]=ttsolve(xc,...
                varargin{:},'mpri',0,'estpol',1);
        end
        
        % align records
        data0=multiply(timeshift(data0,-arr),pol);

        % plot alignment
        set(ax(4),'visible','on');
        info.handles(4)=plot2(normalize(data0),'ax',ax(4));
        
        % force user to decide
        choice=0;
        while(~choice)
            % ask user if they are happy with alignment
            choice=menu('KEEP THIS ALIGNMENT?',...
                'YES',...
                'NO - TRY AGAIN',...
                'NO - TRY AGAIN WITH THESE OFFSETS');
            switch choice
                case 1 % rainbow's end
                    % allow picking the stack to deal with dc-shift
                    [b,e,delta]=getheader(data0,'b','e','delta');
                    set(ax(5),'visible','on');
                    [onset,info.handles(5)]=pickstack(ax(5),...
                        normalize(data0),2/min(delta),[],min(b),max(e),0);
                    data0=timeshift(data0,-onset);
                    arr=-getheader(data0,'o');
                    happy_user=true;
                case 2 % never never quit!
                    data0=data;
                    for a=1:numel(ax)
                        cla(ax(a),'reset');
                    end
                    set(ax,'visible','off','color','k');
                case 3 % standing on the shoulders of those before me
                    info.iteration=info.iteration+1;
                    arr=-getheader(data0,'o');
                    data0=timeshift(data,-o-arr);
                    for a=1:numel(ax)
                        cla(ax(a),'reset');
                    end
                    set(ax,'visible','off','color','k');
            end
        end
    end
    
    % put results in info
    info.solution.arr=arr;
    info.solution.arrerr=err;
    info.solution.pol=pol;
    info.solution.zmean=zmean;
    info.solution.zstd=zstd;
    info.solution.nc=nc;
    
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
