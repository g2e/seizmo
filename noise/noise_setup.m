function []=noise_setup(indir,outdir,varargin)
%NOISE_SETUP    Convert a directory of data to a year/time/file filesystem
%
%    Usage:    noise_setup(indir,outdir)
%              noise_setup(indir,outdir,'opt1',value1,...,'optN',valueN)
%
%    Description:
%     NOISE_SETUP(INDIR,OUTDIR) reads in the records in directory INDIR,
%     merges them where possible and creates a directory system of
%     timesections where each timesection directory contains records 3
%     hours or 180 minutes long starting at the beginning of each day
%     (with no overlap) in OUTDIR.  The directory layout is:
%      OUTDIR/YYYY/YYYY.DAY.HH.MM.SS_YYYY.DAY.HH.MM.SS/RECORDS
%     where YYYY is the year and the two YYYY.DAY.HH.MM.SS strings give the
%     date/time of the timesection start and end respectively.  The A & F
%     fields for each record are set to the timesection time limits and the
%     reference time is set to the timesection start.
%
%     NOISE_SETUP(INDIR,OUTDIR,'OPT1',VALUE1,...,'OPTN',VALUEN) changes the
%     indicated options to another value.  The following options are
%     configurable:
%      LENGTH       - time section length in MINUTES (default is 180)
%      OVERLAP      - time section overlap in MINUTES (default is 0)
%      TIMESTART    - no time sections before this UTC time (default is [])
%      TIMEEND      - no time sections after this UTC time (default is [])
%      LATRNG       - include stations in this latitude range []
%      LONRNG       - include stations in this longitude range []
%      NETWORKS     - include records with these network codes []
%      STATIONS     - include records with these station codes []
%      STREAMS      - include records with these stream codes []
%      COMPONENTS   - include records with these component codes []
%      FILENAMES    - limit processing to files with these filenames []
%      QUIETWRITE   - quietly overwrite OUTDIR (default is false)
%
%     LENGTH & OVERLAP must be whole numbers (eg. 4.5 minutes is NOT
%     allowed!).  Please note that in order to have consistency from day to
%     day the timesection sequence is forced to restart at the start of
%     each day.  This means that timesections near the day boundary may not
%     honor the OVERLAP value.  Also please note that timesections that
%     overlap TIMESTART or TIMEEND are output too.
%
%    Notes:
%     - Handles one component at a time. If all the files for one component
%       of data exceeds your RAM then this will cause trouble. Try working
%       on smaller portions of the input dataset if that is the case.
%     - Please cleanup your headers beforehand!  Make sure that the KNETWK,
%       KSTNM, KHOLE & KCMPNM fields are filled with the appropriate info
%       beforehand (FIX_RDSEED_V48 is helpful)!
%     - The reason that timesections are 3 hours by default is to allow
%       straightforward comparison to the global WW3 model:
%        http://polar.ncep.noaa.gov/waves/viewer.shtml?
%        ftp://polar.ncep.noaa.gov/pub/history/waves/
%        http://polar.ncep.noaa.gov/waves/pres/primer/primer_1.html
%
%    Header changes: NPTS, B, E, DEP*, Z, IZTYPE, A, F
%
%    Examples:
%     % Setup for 25 hour time sections and 1hr overlap:
%     noise_setup('raw','setup-25hour','l',60*25,'o',60)
%
%     % Setup for 15 minute time sections and 80% overlap:
%     noise_setup('some/dir','my_15min','l',15,'o',12)
%
%    See also: NOISE_PROCESS, NOISE_STACK, NOISE_OVERVIEW

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
%        Jan. 13, 2012 - lat/lon, kname, filename subsetting
%        Jan. 23, 2012 - options are now checked, doc update, sync & a/f
%                        marker bugfix
%        Jan. 26, 2012 - use 3hr no-overlap default
%        Jan. 27, 2012 - minor improvement to string option handling
%        Jan. 28, 2012 - checking state bugfix
%        Feb.  1, 2012 - output kname vs progress bar, no iztype warning,
%                        subsetting by user ts/te options, 3char cmp check,
%                        all time operations are in utc
%        Feb.  7, 2012 - merge to meld update
%        Mar.  8, 2012 - drop UTC in timediff where unnecessary
%        Mar. 15, 2012 - parallel verbose fixes
%        Apr. 25, 2012 - minor example fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 25, 2012 at 11:15 GMT

% todo:

% check nargin
error(nargchk(2,inf,nargin));

% parse/check options
opt=noise_setup_parameters(varargin{:});

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
    if(~opt.QUIETWRITE)
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

% parallel processing setup
verbose=seizmoverbose;
%matlabpool(4); % PARALLEL

% read in data headers
if(verbose); disp('READING DATAFILE HEADERS (MAY TAKE A WHILE)!)'); end
data=readheader(strcat(indir,fs,opt.FILENAMES));

% require 3char cmp code
if(size(char(getheader(data,'kcmpnm')),2)~=3)
    error('seizmo:noise_setup:badDATA',...
        'KCMPNM of some records is not 3 characters long!');
end

% limit stations based on user input
if(~isempty(opt.TIMESTART))
    eutc=cell2mat(getheader(data,'e utc'));
    data=data(timediff(opt.TIMESTART,eutc)>0);
    checkdata(data);
end
if(~isempty(opt.TIMEEND))
    butc=cell2mat(getheader(data,'b utc'));
    data=data(timediff(opt.TIMEEND,butc)<0);
    checkdata(data);
end
if(~isempty(opt.LATRNG))
    stla=getheader(data,'stla');
    data=data(stla>=min(opt.LATRNG) & stla<=max(opt.LATRNG));
    checkdata(data);
end
if(~isempty(opt.LONRNG))
    stlo=getheader(data,'stlo');
    data=data(stlo>=min(opt.LONRNG) & stlo<=max(opt.LONRNG));
    checkdata(data);
end
if(~isempty(opt.NETWORKS))
    knetwk=getheader(data,'knetwk');
    data=data(ismember(lower(knetwk),opt.NETWORKS));
    checkdata(data);
end
if(~isempty(opt.STATIONS))
    kstnm=getheader(data,'kstnm');
    data=data(ismember(lower(kstnm),opt.STATIONS));
    checkdata(data);
end
if(~isempty(opt.STREAMS))
    khole=getheader(data,'khole');
    data=data(ismember(lower(khole),opt.STREAMS));
    checkdata(data);
end
if(~isempty(opt.COMPONENTS))
    kcmpnm=getheader(data,'kcmpnm');
    data=data(ismember(lower(kcmpnm),opt.COMPONENTS));
    checkdata(data);
end

% detail message
if(verbose)
    disp('MAKING DIRECTORY STRUCTURE FOR NOISE ANALYSIS');
end

% get name-based indexing
[cmpidx,kname]=getcomponentidx(data);

% loop over components
%parfor i=1:max(cmpidx) % PARALLEL
for i=1:max(cmpidx) % SERIAL
    % force quietness even in parfor (which resets globals)
    seizmoverbose(false);
    
    % detail message
    if(verbose); disp(['PROCESSING: ' kname{i}]); end
    
    % read in component data
    cdata=readdata(data(i==cmpidx));
    
    % merge the component data (does header check)
    checkoperr('invalid_iztype','ignore');
    cdata=meld(cdata);
    checkoperr('invalid_iztype','warn');
        
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
            if((~isempty(opt.TIMESTART) && yr<opt.TIMESTART(1)) ...
                    || (~isempty(opt.TIMEEND) && yr>opt.TIMEEND(1)))
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
                if((~isempty(opt.TIMESTART) ...
                        && datenum(doy2cal([yr day]))...
                        <fix(gregorian2serial(opt.TIMESTART))) ...
                        || (~isempty(opt.TIMEEND) ...
                        && datenum(doy2cal([yr day]))...
                        >fix(gregorian2serial(opt.TIMEEND))))
                    continue;
                end
                
                % loop over time sections
                % 1439 is 1 minute shy of 1 day
                for ts=0:opt.LENGTH-opt.OVERLAP:1439
                    % get time section limits (we avoid utc here)
                    tsbgn=fixtimes([yr day 0 ts 0]);
                    tsend=fixtimes([yr day 0 ts+opt.LENGTH 0]);
                    
                    % skip time section if not in user-defined range
                    if((~isempty(opt.TIMESTART) ...
                            && timediff(opt.TIMESTART,tsend)<=0) ...
                            || (~isempty(opt.TIMEEND) ...
                            && timediff(opt.TIMEEND,tsbgn)>=0))
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
                    tsdata=synchronize(tsdata,tsbgn);
                    tsdata=changeheader(tsdata,'iztype','ia',...
                        'a',0,'f',timediff(tsbgn,tsend,'utc'));
                    
                    % write
                    writeseizmo(tsdata,...
                        'path',[outdir fs num2str(yr) fs tsdir]);
                end
            end
        end
        
        % toggle checking back
        seizmocheck_state(oldseizmocheckstate);
        checkheader_state(oldcheckheaderstate);
    catch
        % toggle checking back
        seizmocheck_state(oldseizmocheckstate);
        checkheader_state(oldcheckheaderstate);
        
        % parallel processing takedown & fix verbosity
        %matlabpool close; % PARALLEL
        seizmoverbose(verbose);
        
        % rethrow error
        error(lasterror);
    end
end

% parallel processing takedown & fix verbosity
%matlabpool close; % PARALLEL
seizmoverbose(verbose);

end


function [opt]=noise_setup_parameters(varargin)
% parses/checks noise_setup pv pairs

% defaults
varargin=[{'l' 180 'o' 0 'q' false 'ts' [] 'te' [] 'lat' [] 'lon' [] ...
    'net' [] 'sta' [] 'str' [] 'cmp' [] 'file' []} varargin];

% require option/value pairs
if(mod(nargin,2))
    error('seizmo:noise_setup:badInput',...
        'Unpaired option/value pair given!');
elseif(~iscellstr(varargin(1:2:end)))
    error('seizmo:noise_setup:badInput',...
        'Options must be specified as strings!');
end

% get user input
for i=1:2:numel(varargin)
    switch lower(varargin{i})
        case {'l' 'len' 'length'}
            if(isempty(varargin{i+1})); continue; end
            opt.LENGTH=varargin{i+1};
        case {'o' 'ov' 'overlap'}
            if(isempty(varargin{i+1})); continue; end
            opt.OVERLAP=varargin{i+1};
        case {'q' 'qw' 'quiet' 'qwrite' 'quietwrite'}
            if(isempty(varargin{i+1})); continue; end
            opt.QUIETWRITE=varargin{i+1};
        case {'ts' 'tstart' 'timestart'}
            opt.TIMESTART=varargin{i+1};
        case {'te' 'tend' 'timeend'}
            opt.TIMEEND=varargin{i+1};
        case {'lat' 'la' 'lar' 'latr' 'larng' 'latitude' 'latrng'}
            opt.LATRNG=varargin{i+1};
        case {'lon' 'lo' 'lor' 'lonr' 'lorng' 'longitude' 'lonrng'}
            opt.LONRNG=varargin{i+1};
        case {'knetwk' 'n' 'net' 'netwk' 'network' 'nets' 'networks'}
            opt.NETWORKS=varargin{i+1};
        case {'kstnm' 'st' 'sta' 'stn' 'stns' 'stations' 'station'}
            opt.STATIONS=varargin{i+1};
        case {'khole' 'hole' 'holes' 'str' 'strs' 'stream' 'streams'}
            opt.STREAMS=varargin{i+1};
        case {'kcmpnm' 'cmpnm' 'cmp' 'cmps' 'component' 'components'}
            opt.COMPONENTS=varargin{i+1};
        case {'f' 'file' 'filename' 'files' 'filenames'}
            opt.FILENAMES=varargin{i+1};
        otherwise
            error('seizmo:noise_setup:badInput',...
                'Unknown Option: %s !',varargin{i});
    end
end

% fix string options to be cellstr vectors
if(ischar(opt.NETWORKS)); opt.NETWORKS=cellstr(opt.NETWORKS); end
if(ischar(opt.STATIONS)); opt.STATIONS=cellstr(opt.STATIONS); end
if(ischar(opt.STREAMS)); opt.STREAMS=cellstr(opt.STREAMS); end
if(ischar(opt.COMPONENTS)); opt.COMPONENTS=cellstr(opt.COMPONENTS); end
if(ischar(opt.FILENAMES)); opt.FILENAMES=cellstr(opt.FILENAMES); end
if(iscellstr(opt.NETWORKS))
    opt.NETWORKS=unique(lower(opt.NETWORKS(:)));
end
if(iscellstr(opt.STATIONS))
    opt.STATIONS=unique(lower(opt.STATIONS(:)));
end
if(iscellstr(opt.STREAMS)); opt.STREAMS=unique(lower(opt.STREAMS(:))); end
if(iscellstr(opt.COMPONENTS))
    opt.COMPONENTS=unique(lower(opt.COMPONENTS(:)));
end
if(iscellstr(opt.FILENAMES)); opt.FILENAMES=unique(opt.FILENAMES(:)); end

% check options
szs=size(opt.TIMESTART);
sze=size(opt.TIMEEND);
if(~isscalar(opt.LENGTH) || ~isreal(opt.LENGTH) || opt.LENGTH<=0)
    error('seizmo:noise_setup:badInput',...
        'LENGTH must be a positive real-valued scalar in minutes!');
elseif(~isscalar(opt.OVERLAP) || ~isreal(opt.OVERLAP)) % negative for gaps
    error('seizmo:noise_setup:badInput',...
        'OVERLAP must be a real-valued scalar in minutes!');
elseif(~isscalar(opt.QUIETWRITE) || ~islogical(opt.QUIETWRITE))
    error('seizmo:noise_setup:badInput',...
        'QUIETWRITE flag must be a scalar logical!');
elseif(~isempty(opt.TIMESTART) && (numel(szs)>2 || szs(1)~=1 ...
        || all(szs(2)~=[2 3 5 6]) || ~isnumeric(opt.TIMESTART) ...
        || ~isreal(opt.TIMESTART)))
    error('seizmo:noise_setup:badInput',...
        'TIMESTART must be a recognized date-time vector!');
elseif(~isempty(opt.TIMEEND) && (numel(sze)>2 || sze(1)~=1 ...
        || all(sze(2)~=[2 3 5 6]) || ~isnumeric(opt.TIMEEND) ...
        || ~isreal(opt.TIMEEND)))
    error('seizmo:noise_setup:badInput',...
        'TIMEEND must be a recognized date-time vector!');
elseif(~isempty(opt.LATRNG) && (~isnumeric(opt.LATRNG) ...
        || ~isreal(opt.LATRNG) || numel(opt.LATRNG)~=2 ...
        || size(opt.LATRNG,2)~=2 || numel(size(opt.LATRNG))~=2))
    error('seizmo:noise_setup:badInput',...
        'LATRNG must be a 2 element numeric vector as [LOW HIGH]!');
elseif(~isempty(opt.LONRNG) && (~isnumeric(opt.LONRNG) ...
        || ~isreal(opt.LONRNG) || numel(opt.LONRNG)~=2 ...
        || size(opt.LONRNG,2)~=2 || numel(size(opt.LONRNG))~=2))
    error('seizmo:noise_setup:badInput',...
        'LONRNG must be a 2 element numeric vector as [LOW HIGH]!');
elseif(~isempty(opt.NETWORKS) && (~iscellstr(opt.NETWORKS)))
    error('seizmo:noise_setup:badInput',...
        'NETWORKS must be a string list of allowed network codes!');
elseif(~isempty(opt.STATIONS) && (~iscellstr(opt.STATIONS)))
    error('seizmo:noise_setup:badInput',...
        'STATIONS must be a string list of allowed station codes!');
elseif(~isempty(opt.STREAMS) && (~iscellstr(opt.STREAMS)))
    error('seizmo:noise_setup:badInput',...
        'STREAMS must be a string list of allowed stream codes!');
elseif(~isempty(opt.COMPONENTS) && (~iscellstr(opt.COMPONENTS)))
    error('seizmo:noise_setup:badInput',...
        'COMPONENTS must be a string list of allowed component codes!');
elseif(~isempty(opt.FILENAMES) && (~iscellstr(opt.FILENAMES)))
    error('seizmo:noise_setup:badInput',...
        'FILENAMES must be a string list of allowed files!');
end

end


function []=checkdata(data)
if(isempty(data))
    error('seizmo:noise_setup:badParam',...
        'No records meet these parameters!');
end
end
