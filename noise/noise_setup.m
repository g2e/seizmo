function []=noise_setup(indir,outdir,varargin)
%NOISE_SETUP    Convert a directory of data to a year/time/file filesystem
%
%    Usage:    noise_setup(indir,outdir)
%              noise_setup(indir,outdir,'opt1',value1,...,'optN',valueN)
%
%    Description:
%     NOISE_SETUP(INDIR,OUTDIR) reads in the records in INDIR, combines
%     them where possible and creates a directory system of "day" records
%     (25 hours long starting at the beginning of each day -- gives 1 hour
%     of overlap) in OUTDIR.  The directory layout is:
%      OUTDIR/YYYY/YYYY.DAY.HH.MM.SS_YYYY.DAY.HH.MM.SS/RECORDS
%     The two times given in the directory structure are the time limits of
%     the time section.
%
%     NOISE_SETUP(INDIR,OUTDIR,'OPT1',VALUE1,...,'OPTN',VALUEN) changes the
%     indicated options to another value.  The following options are
%     configurable:
%      LENGTH       - time section length in MINUTES (default is 1500)
%      OVERLAP      - time section overlap in MINUTES (default is 60)
%      TIMESTART    - no time sections before this time (default is [])
%      TIMEEND      - no time sections after this time (default is [])
%      QUIETWRITE   - quietyly overwrite OUTDIR (default is false)
%     LENGTH & OVERLAP must be whole numbers (eg. 4.5 minutes is NOT
%     allowed!).  In order to have consistency from day to day the time
%     section sequence is forced to restart at the start of each day.  This
%     means that time sections near the day boundary may not honor the
%     OVERLAP value.  Also please note that time sections that overlap
%     TIMESTART or TIMEEND are processed.
%
%    Notes:
%     - Handles one component at a time. If all the files for one component
%       of data exceeds your RAM then this will cause trouble. Try working
%       on smaller portions of the input dataset if that is the case.
%     - Make sure that the KNETWK, KSTNM, KHOLE & KCMPNM fields are filled
%       with the appropriate information beforehand (also probably run
%       FIX_RDSEED_V48/50 or whatever is appropriate)!
%
%    Header changes: NPTS, B, E, DEP*, Z, IZTYPE
%
%    Examples:
%     % Setup for 3 hour time sections and no overlap:
%     noise_setup('raw','setup-3hour','l',60*3,'o',0)
%
%     % Setup for 15 minute time sections and 80% overlap:
%     noise_setup('some/dir','my_15min','l',15,'o',12)
%
%    See also: NOISE_PROCESS, NOISE_STACK, NOISE_WORKFLOW

%     Version History:
%        June 18, 2010 - initial version
%        June 30, 2010 - better catching of errors
%        Sep. 21, 2010 - commented out parallel processing lines
%        Nov.  3, 2011 - allow variable length/overlap timesections (so
%                        this replaces mergecut too), modified timesection
%                        directory naming to include end time, renamed from
%                        DAYDIRS_MAKE to NOISE_SETUP
%        Nov. 22, 2011 - options are now parameter/value pairs, added
%                        timestart & timeend options, updated docs, sync
%                        data to time section start, set a/f markers to
%                        time section limits
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 22, 2011 at 11:15 GMT

% todo:
% - latrng, lonrng, stations

% check nargin
error(nargchk(2,inf,nargin));

% parse/check options
[len,ov,qw,userts,userte]=noise_setup_parameters(varargin{:});

% check directories
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:noise_setup:fileNotString',...
        'INDIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:noise_setup:dirConflict',...
        ['Input Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~ischar(outdir) || ~isvector(outdir))
    error('seizmo:noise_setup:fileNotString',...
        'OUTDIR must be a string!');
end
if(exist(outdir,'file'))
    if(~exist(outdir,'dir'))
        error('seizmo:noise_setup:dirConflict',...
            'Output Directory: %s\nIs a file!',outdir);
    end
    if(~qw)
        fprintf('Output Directory: %s\nDirectory Exists!\n',outdir);
        reply=input('Overwrite? Y/N [N]: ','s');
        if(isempty(reply) || ~strncmpi(reply,'y',1))
            disp('Not overwriting!');
            return;
        end
        disp('Overwriting!');
    end
end

% directory separator
fs=filesep;

% parallel processing setup (up to 8 instances)
%matlabpool(4); % PARALLEL

% verbosity  (turn it off for the loop)
verbose=seizmoverbose(false);

% read in data headers
data=readheader(indir);

% number of records
nrecs=numel(data);

% detail message
if(verbose)
    disp('Making directory structure for noise analysis');
    print_time_left(0,nrecs);
end

% get name-based indexing
[cmpidx,kname]=getcomponentidx(data);
%ncrecs=zeros(max(cmpidx),1); % PARALLEL
%for i=1:max(cmpidx); ncrecs(i)=sum(i==cmpidx); end % PARALLEL

% loop over components
nprocessed=0; % SERIAL
%parfor i=1:max(cmpidx) % PARALLEL
for i=1:max(cmpidx) % SERIAL
    % read in component data
    cdata=readdata(data(i==cmpidx));
    ncrecs=numel(cdata); % SERIAL
    
    % merge the component data (does header check)
    cdata=merge(cdata);
        
    % turn off checking
    oldseizmocheckstate=seizmocheck_state(false);
    oldcheckheaderstate=checkheader_state(false);
    
    try
        % get absolute time limits
        [butc,eutc]=getheader(cdata,'b utc','e utc');
        butc=cell2mat(butc); eutc=cell2mat(eutc);
        minyr=min(butc(:,1)); maxyr=max(eutc(:,1));
        minday=min(butc(butc(:,1)==minyr,2));
        maxday=max(eutc(eutc(:,1)==maxyr,2));
        
        % loop over years
        for yr=minyr:maxyr
            % skip year if not in user-defined range
            if((~isempty(userts) && yr<userts(1)) ...
                    || (~isempty(userte) && yr>userte(1)))
                continue;
            end
            
            % which days?
            if(minyr==maxyr)
                days=minday:maxday;
            elseif(yr==minyr)
                days=minday:365+isleapyear(minyr);
            elseif(yr==maxyr)
                days=1:maxday;
            else
                days=1:365+isleapyear(yr);
            end
            
            % loop over the days
            for day=days
                % skip day if not in user-defined range
                if((~isempty(userts) ...
                        && datenum(doy2cal([yr day]))...
                        <fix(gregorian2serial(userts))) ...
                        || (~isempty(userte) ...
                        && datenum(doy2cal([yr day]))...
                        >fix(gregorian2serial(userte))))
                    continue;
                end
                
                % loop over time sections
                % 1439 is 1 minute shy of 1 day
                for ts=0:len-ov:1439
                    % get time section limits (we avoid utc here)
                    tsbgn=fixtimes([yr day 0 ts 0]);
                    tsend=fixtimes([yr day 0 ts+len 0]);
                    
                    % skip time section if not in user-defined range
                    if((~isempty(userts) ...
                            && timediff(userts,tsend)<=0) ...
                            || (~isempty(userte) ...
                            && timediff(userte,tsbgn)>=0))
                        continue;
                    end
                    
                    % time section directory name
                    tsdir=sprintf(['%04d.%03d.%02d.%02d.%02d_' ...
                        '%04d.%03d.%02d.%02d.%02d'],tsbgn,tsend);
                    
                    % debugging
                    %disp([kname{i} ' ' num2str(yr) ' ' tsdir])
                    
                    % skip if time section contains no records
                    if(all(timediff(butc,tsend)<=0 ...
                            | timediff(eutc,tsbgn)>=0))
                        continue;
                    end
                    
                    % cut
                    tsdata=cut(cdata,tsbgn,tsend);
                    
                    % skip if none
                    if(isempty(tsdata)); continue; end
                    
                    % sync & insert markers for time section
                    % - note we allow some minutes to have leapseconds
                    data=synchronize(data,tsbgn);
                    data=changeheader(data,'iztype','ia',...
                        'a',0,'f',timediff(tsbgn,tsend,'utc'));
                    
                    % write
                    writeseizmo(tsdata,...
                        'path',[outdir fs num2str(yr) fs tsdir]);
                end
            end
        end
    catch
        % toggle checking back
        seizmocheck_state(oldseizmocheckstate);
        checkheader_state(oldcheckheaderstate);
        seizmoverbose(verbose);
        
        % rethrow error
        error(lasterror);
    end
    
    % detail message
    nprocessed=nprocessed+ncrecs; % SERIAL
    if(verbose); print_time_left(nprocessed,nrecs); end % SERIAL
    %if(verbose); print_time_left(sum(ncrecs(1:i)),nrecs); end % PARALLEL
end

% parallel processing takedown & fix verbosity
%matlabpool close; % PARALLEL
seizmoverbose(verbose);

end


function [len,ov,qw,userts,userte]=noise_setup_parameters(varargin)

% defaults
varargin=[{'l' 1500 'o' 60 'q' false 'ts' [] 'te' []} varargin];

% require option/value pairs
if(mod(nargin,2))
    error('seizmo:noise_setup:badInput',...
        'Unpaired option/value pair given!');
elseif(~iscellstr(varargin(1:2:end)))
    error('seizmo:noise_setup:badInput',...
        'Options must be specified as strings!');
end

% get user input
for i=1:2:nargin
    switch lower(varargin{i})
        case {'l' 'len' 'length'}
            if(isempty(varargin{i+1})); continue; end
            len=varargin{i+1};
        case {'o' 'ov' 'overlap'}
            if(isempty(varargin{i+1})); continue; end
            ov=varargin{i+1};
        case {'q' 'qw' 'quiet' 'qwrite' 'quietwrite'}
            if(isempty(varargin{i+1})); continue; end
            qw=varargin{i+1};
        case {'ts' 'tstart' 'timestart'}
            userts=varargin{i+1};
        case {'te' 'tend' 'timeend'}
            userte=varargin{i+1};
        otherwise
            error('seizmo:noise_setup:badInput',...
                'Unknown Option: %s !',varargin{i});
    end
end

% check options
szs=size(userts);
sze=size(userte);
if(~isscalar(len) || ~isreal(len) || len<=0)
    error('seizmo:noise_setup:badInput',...
        'LENGTH must be a positive real-valued scalar in minutes!');
elseif(~isscalar(ov) || ~isreal(ov)) % negative for gaps
    error('seizmo:noise_setup:badInput',...
        'OVERLAP must be a real-valued scalar in minutes!');
elseif(~isscalar(qw) || ~islogical(qw))
    error('seizmo:noise_setup:badInput',...
        'QUIETWRITE flag must be a scalar logical!');
elseif(numel(szs)>2 || szs(1)~=1 || all(szs(2)~=[2 3 5 6]) ...
        || ~isnumeric(userts) || ~isreal(userts))
    error('seizmo:noise_setup:badInput',...
        'TIMESTART must be a recognized date-time vector!');
elseif(numel(sze)>2 || sze(1)~=1 || all(sze(2)~=[2 3 5 6]) ...
        || ~isnumeric(userte) || ~isreal(userte))
    error('seizmo:noise_setup:badInput',...
        'TIMEEND must be a recognized date-time vector!');
end

end
