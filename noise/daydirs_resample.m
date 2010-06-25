function []=daydirs_resample(indir,outdir,rate,o)
%DAYDIRS_RESAMPLE    Sample records in day directories to a new samplerate
%
%    Usage:    daydirs_resample(indir,outdir)
%              daydirs_resample(indir,outdir,samplerate)
%              daydirs_resample(indir,outdir,samplerate,overwrite)
%
%    Description: DAYDIRS_RESAMPLE(INDIR,OUTDIR) resamples records in a day
%     directory layout under directory INDIR to 1 sample per second.  The
%     resampled records are output in an equivalent directory layout under
%     OUTDIR.
%
%     DAYDIRS_RESAMPLE(INDIR,OUTDIR,SR) resamples the samplerate SR.  The
%     default value is 1.
%
%     DAYDIRS_RESAMPLE(INDIR,OUTDIR,SR,OVERWRITE) quietly overwrites
%     pre-existing records in OUTDIR when OVERWRITE is set to TRUE.  By
%     default OVERWRITE is FALSE.
%
%    Notes:
%
%    Header changes: E, NPTS, DELTA, DEP*
%
%    Examples:
%
%    See also: DAYDIRS_MERGECUT_25HRS, DAYDIRS_RINST, DAYDIRS_NORMALIZE,
%              DAYDIRS_CORRELATE, DAYDIRS_ROTCORR, DAYDIRS_STACKCORR,
%              DAYDIRS_MAKE

%     Version History:
%        June 18, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 18, 2010 at 12:55 GMT

% todo:

% check nargin
error(nargchk(2,4,nargin));

% defaults
if(nargin<3 || isempty(rate)); rate=1; end
if(nargin<4 || isempty(o)); o=false; end
if(~isscalar(rate) || ~isreal(rate) || rate<=0)
    error('seizmo:daydirs_resample:badInput',...
        'RATE must be a positive real-valued scalar!');
end
if(~isscalar(o) || ~islogical(o))
    error('seizmo:daydirs_resample:badInput',...
        'OVERWRITE flag must be a scalar logical!');
end

% check directories
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:daydirs_resample:fileNotString',...
        'INDIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:daydirs_resample:dirConflict',...
        ['Input Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~ischar(outdir) || ~isvector(outdir))
    error('seizmo:daydirs_resample:fileNotString',...
        'OUTDIR must be a string!');
end
if(exist(outdir,'file'))
    if(~exist(outdir,'dir'))
        error('seizmo:daydirs_resample:dirConflict',...
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

% get directory separator
fs=filesep;

% parallel processing setup (8 instances)
matlabpool(8);

% get year directories and day directories
dirs=xdir([indir fs]);
dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dirs
years=str2double({dirs.name});
nyears=numel(years);
if(any(isnan(years)))
    error('seizmo:daydirs_resample:badLayout',...
        'Improper directory layout!');
end
jdays=cell(size(years));
for i=1:nyears
    % get day directories
    dirs=xdir([indir fs num2str(years(i))]);
    dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dir
    jdays{i}=str2double({dirs.name});
    if(any(isnan(jdays{i})))
        error('seizmo:daydirs_resample:badLayout',...
            'Improper directory layout!');
    end
end

% verbosity (turn it off for the loop)
verbose=seizmoverbose(false);
if(verbose); disp('Altering samplerate of record(s)'); end

% loop over years
for i=1:nyears
    % working year
    yr=years(i);
    syr=num2str(yr);
    
    % loop over days
    parfor j=1:numel(jdays{i})
        % working julian day
        jday=jdays{i}(j);
        sjday=num2str(jday,'%03d');
        
        % detail message
        if(verbose); disp(['PROCESSING DAY ' syr '.' sjday]); end
        
        % attempt resample
        try
            % read in data
            try
                data=readseizmo([indir fs syr fs sjday fs]);
            catch
                % empty day
                continue;
            end
            
            % resample
            data=syncrates(data,rate);
            
            % write
            writeseizmo(data,'pathchange',{indir outdir});
        catch
            % close pool & fix verbosity
            matlabpool close;
            seizmoverbose(verbose);
            
            % ???
            error(lasterror);
        end
    end
end

% parallel processing takedown
matlabpool close;
seizmoverbose(verbose);

end
