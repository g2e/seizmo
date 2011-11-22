function []=noise_setup(indir,outdir,length,overlap,o)
%NOISE_SETUP    Convert a directory of data to a year/time/file filesystem
%
%    Usage:    noise_setup(indir,outdir)
%              noise_setup(indir,outdir,length,overlap)
%              noise_setup(indir,outdir,length,overlap,overwrite)
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
%     NOISE_SETUP(INDIR,OUTDIR,LENGTH,OVERLAP) allows redefining the time
%     section length LENGTH and overlap between time sections OVERLAP.
%     Both are in __minutes__!  The default LENGTH is 1500 (25 hours) and
%     the default OVERLAP is 60 (1 hour). Non-integer values are not
%     allowed!  In order to have consistency from day to day the time
%     sections are forced to restart at the start of each day.  This means
%     that time sections near the day boundary may not honor the OVERLAP
%     value.
%
%     NOISE_SETUP(INDIR,OUTDIR,LENGTH,OVERLAP,OVERWRITE) quietly writes
%     in OUTDIR when OVERWRITE is set to TRUE.  By default OVERWRITE is
%     FALSE.
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
%     noise_setup('raw','setup-3hour',60*3,0)
%
%     % Setup for 15 minute time sections and 80% overlap:
%     noise_setup('some/dir','my_15min',15,12)
%
%    See also: NOISE_PROCESS, NOISE_STACK

%     Version History:
%        June 18, 2010 - initial version
%        June 30, 2010 - better catching of errors
%        Sep. 21, 2010 - commented out parallel processing lines
%        Nov.  3, 2011 - allow variable length/overlap timesections (so
%                        this replaces mergecut too), modified timesection
%                        directory naming to include end time, renamed from
%                        DAYDIRS_MAKE to NOISE_SETUP
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  3, 2011 at 11:15 GMT

% todo:
% - time window option

% check nargin
error(nargchk(2,5,nargin));

% defaults
if(nargin<3 || isempty(length)); length=1500; end
if(nargin<4 || isempty(overlap)); overlap=60; end
if(nargin<5 || isempty(o)); o=false; end

% check options
if(~isscalar(length) || ~isreal(length) || length<=0)
    error('seizmo:noise_setup:badInput',...
        'LENGTH must be a positive real-valued scalar in minutes!');
elseif(~isscalar(overlap) || ~isreal(length)) % negative for gaps
    error('seizmo:noise_setup:badInput',...
        'OVERLAP must be a real-valued scalar in minutes!');
elseif(~isscalar(o) || ~islogical(o))
    error('seizmo:noise_setup:badInput',...
        'OVERWRITE flag must be a scalar logical!');
end

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
    if(~o)
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
                % loop over time sections
                % 1439 is 1 minute shy of 1 day
                for ts=0:length-overlap:1439
                    % get time section limits (we avoid utc here)
                    tsbgn=fixtimes([yr day 0 ts 0]);
                    tsend=fixtimes([yr day 0 ts+length 0]);
                    
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
