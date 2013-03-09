function [info,xc,data0]=useralign_quiet(data,varargin)
%USERALIGN_QUIET    NONinteractive alignment of a signal for SEIZMO records
%
%    Usage:    [info,xc,data0]=useralign_quiet(data)
%              [...]=useralign_quiet(data,'option1',value1,...)
%
%    Description:
%     [INFO,XC,DATA0]=USERALIGN_QUIET(DATA) aligns records in SEIZMO struct
%     DATA on a particular signal noninteractively.  This provides an
%     automated alternative to USERALIGN.  The output INFO is setup to
%     appear just like from USERALIGN even though interactive functions
%     were not run (USERMOVEOUT, USERWINDOW, USERTAPER are included and
%     contain the default info).  The final solution is within
%     INFO.SOLUTION.  XC is the struct from CORRELATE reordered by TTSOLVE
%     (contains the correlogram peaks used to find the solution).  DATA0
%     is the processed dataset (windowed, tapered & aligned).  See the next
%     usage form to supply parameters to control windowing, tapering,
%     correlation, and alignment.
%
%     [...]=USERALIGN_QUIET(DATA,'OPTION1',VALUE1,...) allows passing
%     options on to CUT, TAPER, CORRELATE & TTSOLVE for enhanced control of
%     the workflow.  Currently available are all options in TTSOLVE, 3
%     options in CORRELATE ('LAGS', 'NPEAKS', & 'ABSXC'), as well as
%     'WINDOW' & 'TAPER'.  'WINDOW' specifies the window passed to CUT and
%     should be a 1x2 vector of [START END].  'TAPER' specifies the taper
%     width passed to TAPER.  See those functions for details about the
%     options.  Note that NPEAKS must be >=1!
%
%    Notes:
%
%    Examples:
%     % I like to give a SNR estimate of the waveforms to provide more info
%     % for the solver to decide which measurements are more reliable:
%     data=timeshift(data,-findpicks(data,'Pdiff',1));
%     snr=quicksnr(data,[-100 -10],[-10 60]);
%     [info,xc,data0]=useralign_quiet(data,'snr',snr);
%
%    See also: USERALIGN, TTSOLVE, CORRELATE, USERWINDOW, USERTAPER,
%              USERRAISE, USERMOVEOUT, MULTIBANDALIGN

%     Version History:
%        Jan. 14, 2011 - initial version
%        Apr. 16, 2011 - add lags option for correlate
%        Mar. 16, 2012 - doc update
%        Jan. 30, 2013 - update for new correlate (no spacing option now)
%        Feb. 27, 2013 - bugfix: needs 'reltime' correlate option, added
%                        info.solution.origarr for debugging purposes
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 27, 2013 at 11:00 GMT

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
        error('seizmo:useralign_quiet:badInput',...
            'Unpaired OPTION/VALUE!');
    elseif(~iscellstr(varargin(1:2:end)))
        error('seizmo:useralign_quiet:badInput',...
            'All OPTIONs must be specified with a string!');
    end
    
    % default options
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
                    error('seizmo:useralign_quiet:badInput',...
                        'LAGS must be given as [MINLAG MAXLAG]!');
                end
                info.correlate.lags=varargin{i+1};
                keep(i:i+1)=false;
            case 'npeaks'
                if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}) ...
                        || varargin{i+1}~=fix(varargin{i+1}) ...
                        || varargin{i+1}<1)
                    error('seizmo:useralign_quiet:badInput',...
                        'NPEAKS must be an integer >=1 !');
                end
                info.correlate.npeaks=varargin{i+1};
                keep(i:i+1)=false;
            case 'absxc'
                if(~isscalar(varargin{i+1}) || ~islogical(varargin{i+1}))
                    error('seizmo:useralign_quiet:badInput',...
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
    outc=cell2mat(getheader(data,'o utc'));
    if(any(abs(timediff(outc(1,:),outc))>0.002))
        error('seizmo:useralign_quiet:oFieldVaries',...
            'O field must correspond to one UTC time for all records!');
    end
    
    % copy data to working dataset (so user can go back)
    data0=data;
    
    % some fields that won't change
    info.iteration=1;
    info.handles=-1*ones(1,5);
    
    % "usermoveout"
    info.usermoveout.moveout=0;
    info.usermoveout.adjust=zeros(numel(data0),1);
    
    % "userwindow"
    info.userwindow.func=@deal;
    info.userwindow.fill=false;
    data0=cut(data0,...
        info.userwindow.limits(1),info.userwindow.limits(2));
    
    % "usertaper"
    info.usertaper.type=[];
    info.usertaper.option=[];
    info.usertaper.func=@deal;
    data0=taper(data0,info.usertaper.width);
    
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
    
    % debugging
    %figure;
    %ax=makesubplots(3,1,1:3);
    %plot2(normalize(data0),'ax',ax(1),'bgc','w');
    %plot2(normalize(timeshift(data0,-arr)),'ax',ax(2),'bgc','w');
    %plot2(normalize(timeshift(data0,arr)),'ax',ax(3),'bgc','w');
    info.solution.origarr=arr;
    
    % align records
    data0=multiply(timeshift(data0,-arr),pol);
    arr=-getheader(data0,'o');
    
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
