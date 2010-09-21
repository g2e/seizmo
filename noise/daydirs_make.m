function []=daydirs_make(indir,outdir,o)
%DAYDIRS_MAKE    Convert a directory of data to a dir/year/day/file system
%
%    Usage:    daydirs_make(indir,outdir)
%              daydirs_make(indir,outdir,overwrite)
%
%    Description: DAYDIRS_MAKE(INDIR,OUTDIR) takes SEIZMO records in
%     directory INDIR and cuts them up into day-long records in a directory
%     layout under OUTDIR.  The directory layout goes:
%      OUTDIR/YEAR/JULDAY/RECORDS
%
%     DAYDIRS_MAKE(INDIR,OUTDIR,OVERWRITE) quietly overwrites pre-existing
%     records in OUTDIR when OVERWRITE is set to TRUE.  By default
%     OVERWRITE is FALSE.
%
%    Notes:
%     - Not setup to handle unevenly sampled records!
%     - Don't forget to run FIX_RDSEED_V48/50 prior to this!
%
%    Header changes: NPTS, B, E, DEP*, Z (Set to beginning of day), IZTYPE
%
%    Examples:
%
%    See also: DAYDIRS_MERGECUT_25HRS, DAYDIRS_RINST, DAYDIRS_NORMALIZE,
%              DAYDIRS_CORRELATE, DAYDIRS_ROTCORR, DAYDIRS_STACKCORR,
%              DAYDIRS_RESAMPLE

%     Version History:
%        June 18, 2010 - initial version
%        June 30, 2010 - better catching of errors
%        Sep. 21, 2010 - commented out parallel processing lines
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 21, 2010 at 11:15 GMT

% todo:
% - there appears to be a leapsecond problem in this that I haven't solved

% check nargin
error(nargchk(2,3,nargin));

% default overwrite to false
if(nargin==2 || isempty(o)); o=false; end
if(~isscalar(o) || ~islogical(o))
    error('seizmo:daydirs_make:badInput',...
        'OVERWRITE flag must be a scalar logical!');
end

% check directories
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:daydirs_make:fileNotString',...
        'INDIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:daydirs_make:dirConflict',...
        ['Input Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~ischar(outdir) || ~isvector(outdir))
    error('seizmo:daydirs_make:fileNotString',...
        'OUTDIR must be a string!');
end
if(exist(outdir,'file'))
    if(~exist(outdir,'dir'))
        error('seizmo:daydirs_make:dirConflict',...
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

% parallel processing setup (8 instances)
%matlabpool(8);

% read in data headers
data=readheader(indir);

% number of records
nrecs=numel(data);

% get absolute time limits (& reference time)
[butc,eutc,delta,z]=getheader(data,'b utc','e utc','delta','z');
butc=cell2mat(butc);
eutc=cell2mat(eutc);
z=cell2mat(z);

% change reference time type
data=changeheader(data,'iztype','iunkn');

% verbosity (turn it off for the loop)
verbose=seizmoverbose(false);
if(verbose)
    disp('Making directory structure of day records');
    print_time_left(0,nrecs);
end

% loop over every record
for i=1:nrecs
%parfor i=1:nrecs
    % read in record
    rec=readdata(data(i));
    
    % get days associated with this record
    days=fix(gregorian2serial([butc(i,:); eutc(i,:)]));
    days=days(1):days(2);
    ndays=numel(days);
    dates=serial2gregorian(days,'doytime');
    
    % day window limits (relative to record's reference time)
    % - from first possible sample of the day to last
    ts=timediff(z(i,:),dates,'utc');
    te=timediff(z(i,:),dates+repmat([0 1 0 0 -delta(i)],ndays,1),'utc');
    
    % loop over days 
    for j=1:ndays
        try
            % cut out day
            dayrec=cut(rec,ts(j),te(j));
            
            % set reference time to start of day
            dayrec=timeshift(dayrec,-ts(j),'iday');
            
            % write
            writeseizmo(dayrec,'path',...
                [outdir fs num2str(dates(j,1)) fs ...
                num2str(dates(j,2),'%03d') fs]);
        catch
            disp(i)
            disp(z(i,:))
            disp(data(i).name)
            disp(j)
            disp(days(j))
            disp(dates(j,:))
            disp(ts(j))
            disp(te(j))
            bad=lasterror;
            warning(bad.message)    
        end
    end
    
    % detail message
    if(verbose); print_time_left(i,nrecs); end
end

% parallel processing takedown & fix verbosity
%matlabpool close;
seizmoverbose(verbose);

end
